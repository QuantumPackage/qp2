#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/alloc.h>
#include <caml/custom.h>
#include <caml/threads.h>

#include <string.h>



/* Adapted from
   https://github.com/monadbobo/ocaml-core/blob/master/base/core/lib/linux_ext_stubs.c
*/

#include <unistd.h>
#include <sys/ioctl.h>
#include <net/if.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

CAMLprim value get_ipv4_address_for_interface(value v_interface)
{
  CAMLparam1(v_interface);
  struct ifreq ifr;
  int fd = -1;
  value res;
  char* error = NULL;

  memset(&ifr, 0, sizeof(ifr));
  ifr.ifr_addr.sa_family = AF_INET;
  /* [ifr] is already initialized to zero, so it doesn't matter if the
     incoming string is too long, and [strncpy] fails to add a \0. */
  strncpy(ifr.ifr_name, String_val(v_interface), IFNAMSIZ - 1);

  caml_enter_blocking_section();
  fd = socket(AF_INET, SOCK_DGRAM, 0);

  if (fd == -1)
    error = "error: couldn't allocate socket";
  else {
    if (ioctl(fd, SIOCGIFADDR, &ifr) < 0)
      error = "error: ioctl(fd, SIOCGIFADDR, ...) failed";

    (void) close(fd);
  }

  caml_leave_blocking_section();

  if (error == NULL) {
    /* This is weird but doing the usual casting causes errors when using
     * the new gcc on CentOS 6.  This solution was picked up on Red Hat's
     * bugzilla or something.  It also works to memcpy a sockaddr into
     * a sockaddr_in.  This is faster hopefully.
     */
    union {
      struct sockaddr sa;
      struct sockaddr_in sain;
    } u;
    u.sa = ifr.ifr_addr;
    res = caml_copy_string(inet_ntoa(u.sain.sin_addr));
  }
  else
    res = caml_copy_string(error);
  CAMLreturn(res);

}
 

