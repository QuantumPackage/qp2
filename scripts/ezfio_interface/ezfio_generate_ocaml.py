#!/usr/bin/env python2
"""
This program generates all the
OCaml templates needed by qp_edit

You can see `ezfio_generate_provider.py`
for an example.
"""

import sys

# If type in **kwargs
from ei_handler import Type


class EZFIO_ocaml(object):

    def __init__(self, **kwargs):

        for k, v in kwargs.iteritems():

            try:
                exec "self.{0} = {1}".format(k, v)
            except NameError:
                exec "self.{0} = '{1}'".format(k, v)

    @property
    def Ocaml_type(self):
        return self.type.ocaml.capitalize()

    @property
    def ocaml_type(self):
        return self.type.ocaml

    @property
    def fancy_type(self):
        return self.type.fancy

    def check_if_init(self, l_arg, name):
        for i in l_arg:
            try:
                exec "self.{0}".format(i)
            except AttributeError:
                msg = "You need to provide a '{0}' for creating {1}"
                raise KeyError(msg.format(i, name))

    def create_read(self):
        '''
        You need to instantiate the EZFIO_ocaml with this keyword argument
        ezfio_dir  = str
        ezfio_name = str
        type       = Named_tuple(fancy_type, ocaml_type, fortrant_type)

        Return the read template
        '''

        # ~#~#~#~#~#~#~#~ #
        # C h e c k i n g #
        # ~#~#~#~#~#~#~#~ #

        self.check_if_init(["ezfio_dir", "ezfio_name", "type"],
                           sys._getframe().f_code.co_name)

        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
        # C r e a t e _ t e m pl a t e #
        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #

        l_template = [
            '(* Read snippet for {self.ezfio_name} *)',
            'let read_{self.ezfio_name} () =',
            '  if not (Ezfio.has_{self.ezfio_dir}_{self.ezfio_name} ()) then',
            '     get_default "{self.ezfio_name}"',
            '     |> {self.Ocaml_type}.of_string',
            '     |> Ezfio.set_{self.ezfio_dir}_{self.ezfio_name}',
            '  ;',
            '  Ezfio.get_{self.ezfio_dir}_{self.ezfio_name} ()']

        if self.fancy_type:
            l_template += ["    |> {self.fancy_type}.of_{self.ocaml_type}"]

        l_template += [";;"]

        template = "\n  ".join(l_template)

        # ~#~#~#~#~#~ #
        # R e n d e r #
        # ~#~#~#~#~#~ #

        template_rendered = template.format(**locals())

        # ~#~#~#~#~#~ #
        # R e t u r n #
        # ~#~#~#~#~#~ #
        return template_rendered

    def create_write(self):
        '''
        You need to instantiate the EZFIO_ocaml with this keyword argument
        ezfio_dir  = str
        ezfio_name = str
        type       = Named_tuple(fancy_type, ocaml_type, fortrant_type)

        Return the read template
        '''

        # ~#~#~#~#~#~#~#~ #
        # C h e c k i n g #
        # ~#~#~#~#~#~#~#~ #

        self.check_if_init(["ezfio_dir", "ezfio_name", "type"],
                           sys._getframe().f_code.co_name)

        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
        # C r e a t e _ t e m pl a t e #
        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #

        l_template = ['(* Write snippet for {self.ezfio_name} *)']

        if self.fancy_type:
            l_template += ['let write_{self.ezfio_name} var = ',
                           '  {self.fancy_type}.to_{self.ocaml_type} var',
                           '  |> Ezfio.set_{self.ezfio_dir}_{self.ezfio_name}']
        else:
            l_template += ['let write_{self.ezfio_name} =',
                           '   Ezfio.set_{self.ezfio_dir}_{self.ezfio_name}']

        l_template += [';;']

        template = "\n  ".join(l_template)

        # ~#~#~#~#~#~ #
        # R e n d e r #
        # ~#~#~#~#~#~ #

        template_rendered = template.format(**locals())

        # ~#~#~#~#~#~ #
        # R e t u r n #
        # ~#~#~#~#~#~ #
        return template_rendered

    def create_type(self):
        '''
        You need to instantiate the EZFIO_ocaml with this keyword argument
        l_ezfio_name  = [provider_name, ...]
        l_type =  [Named_tuple(fancy_type, ocaml_type, fortrant_type), ...]

        Return the type template
        '''

        # ~#~#~#~#~#~#~#~ #
        # C h e c k i n g #
        # ~#~#~#~#~#~#~#~ #

        self.check_if_init(["l_ezfio_name", "l_type"],
                           sys._getframe().f_code.co_name)

        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
        # C r e a t e _ t e m pl a t e #
        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #

        l_template = ["(* Generate type *)",
                      "type t = ",
                      "  {"]

        for p, t in zip(self.l_ezfio_name, self.l_type):

            if t.fancy:
                l_template += ["    {0:<30} : {1}.t;".format(p, t.fancy)]
            else:
                l_template += ["    {0:<30} : {1};".format(p, t.ocaml)]

        l_template += ["  } [@@deriving sexp]",
                       ";;"]

        # ~#~#~#~#~#~ #
        # R e t u r n #
        # ~#~#~#~#~#~ #
        return "\n   ".join(l_template)

    def create_read_global(self):
        '''
        You need to instantiate the EZFIO_ocaml with this keyword argument
        l_ezfio_name  = [ezfio_name, ...]

        Return the read_global template
        '''
        # ~#~#~#~#~#~#~#~ #
        # C h e c k i n g #
        # ~#~#~#~#~#~#~#~ #

        self.check_if_init(["l_ezfio_name"],
                           sys._getframe().f_code.co_name)

        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
        # C r e a t e _ t e m pl a t e #
        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #

        l_template = ["(* Read all *)",
                      "let read() = ",
                      "  Some",
                      "  {"]

        l_template += ["    {0:<30} = read_{0} ();".format(p)
                       for p in self.l_ezfio_name]

        l_template += ["  }",
                       ";;"]

        # ~#~#~#~#~#~ #
        # R e t u r n #
        # ~#~#~#~#~#~ #
        return "\n   ".join(l_template)

    def create_write_global(self):
        '''
        You need to instantiate the EZFIO_ocaml with this keyword argument
        l_ezfio_name  = [provider_name, ...]

        Return the type template
        '''
        # ~#~#~#~#~#~#~#~ #
        # C h e c k i n g #
        # ~#~#~#~#~#~#~#~ #

        self.check_if_init(["l_ezfio_name"],
                           sys._getframe().f_code.co_name)

        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
        # C r e a t e _ t e m pl a t e #
        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #

        l_template = ["(* Write all *)",
                      "let write{ "]
        l_template += ["           {0};".format(p) for p in self.l_ezfio_name]
        l_template += ["         } ="]
        l_template += ["  write_{0:<30} {0};".format(p)
                       for p in self.l_ezfio_name]
        l_template += [";;"]

        # ~#~#~#~#~#~ #
        # R e t u r n #
        # ~#~#~#~#~#~ #
        return "\n   ".join(l_template)

    def create_to_string(self):
        '''
        You need to instantiate the EZFIO_ocaml with this keyword argument
        l_ezfio_name  = [provider_name, ...]
        l_type =  [Named_tuple(fancy_type, ocaml_type, fortrant_type), ...]

        Return the type template
        '''
        # ~#~#~#~#~#~#~#~ #
        # C h e c k i n g #
        # ~#~#~#~#~#~#~#~ #

        self.check_if_init(["l_ezfio_name", "l_type"],
                           sys._getframe().f_code.co_name)

        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
        # C r e a t e _ t e m pl a t e #
        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #

        l_template = ['(* to_string*)',
                      'let to_string b =',
                      '  Printf.sprintf "']

        l_template += ["{0} = %s".format(p) for p in self.l_ezfio_name]
        l_template += ['"']

        for p, t in zip(self.l_ezfio_name, self.l_type):

            if t.fancy:
                str_ = t.fancy
            else:
                str_ = t.ocaml.capitalize()

            l_template += ["    ({0}.to_string b.{1})".format(str_, p)]

        l_template += [";;"]

        # ~#~#~#~#~#~ #
        # R e t u r n #
        # ~#~#~#~#~#~ #
        return "\n   ".join(l_template)

    def create_to_rst(self):
        '''
        You need to instantiate the EZFIO_ocaml with this keyword argument
        l_ezfio_name  = [provider_name, ...]
        l_type =  [Named_tuple(fancy_type, ocaml_type, fortrant_type), ...]

        Return the type template
        '''
        # ~#~#~#~#~#~#~#~ #
        # C h e c k i n g #
        # ~#~#~#~#~#~#~#~ #

        self.check_if_init(["l_ezfio_name", "l_type", "l_doc"],
                           sys._getframe().f_code.co_name)

        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
        # C r e a t e _ t e m pl a t e #
        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #

        l_template = ['(* to_rst*)',
                      'let to_rst b =',
                      '  Printf.sprintf "']

        for p, d in zip(self.l_ezfio_name, self.l_doc):

            l_template += ["{0} ::".format(d),
                           "",
                           "  {0} = %s".format(p),
                           ""]
        l_template += ['"']

        for p, t in zip(self.l_ezfio_name, self.l_type):

            if t.fancy:
                str_ = t.fancy
            else:
                str_ = t.ocaml.capitalize()

            l_template += ["    ({0}.to_string b.{1})".format(str_, p)]

        l_template += ["|> Rst_string.of_string",
                       ";;"]

        # ~#~#~#~#~#~ #
        # R e t u r n #
        # ~#~#~#~#~#~ #
        return "\n   ".join(l_template)

    def create_input_auto_generated(self):
        """
        Generate the include of all the Input_module.lower template
        """

        # ~#~#~#~#~#~#~#~ #
        # C h e c k i n g #
        # ~#~#~#~#~#~#~#~ #

        self.check_if_init(["l_module_lower"],
                           sys._getframe().f_code.co_name)

        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
        # C r e a t e _ t e m pl a t e #
        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
        l_template = ["(* Auto generate. Do not modify by Hand *)"]

        l_template = ['open Qputils;;',
                      'open Qptypes;;',
                      'open Core;;',
                      '']

        for m in self.l_module_lower:
            l_template += ["include Input_{0}".format(m)]

        # ~#~#~#~#~#~ #
        # R e t u r n #
        # ~#~#~#~#~#~ #

        return "\n".join(l_template)

    def create_qp_keywords(self):
        """
        Generate keywords template
        """

        # ~#~#~#~#~#~#~#~ #
        # C h e c k i n g #
        # ~#~#~#~#~#~#~#~ #

        self.check_if_init(["l_module_lower"],
                           sys._getframe().f_code.co_name)

        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
        # C r e a t e _ t e m pl a t e #
        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #

        l_template = ["| {0}".format(m.capitalize())
                      for m in self.l_module_lower]

        # ~#~#~#~#~#~ #
        # R e t u r n #
        # ~#~#~#~#~#~ #

        return "\n".join(l_template)

    def create_qp_keywords_to_string(self):
        """
        Generate keywords to string template
        """

        # ~#~#~#~#~#~#~#~ #
        # C h e c k i n g #
        # ~#~#~#~#~#~#~#~ #

        self.check_if_init(["l_module_lower"],
                           sys._getframe().f_code.co_name)

        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
        # C r e a t e _ t e m pl a t e #
        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #

        l_template = ['| {0} -> "{0}"'.format(m.capitalize())
                      for m in self.l_module_lower]

        # ~#~#~#~#~#~ #
        # R e t u r n #
        # ~#~#~#~#~#~ #

        return "\n".join(l_template)

    def create_qp_section_to_rst(self):
        """
        Generate section to rst
        """
        # ~#~#~#~#~#~#~#~ #
        # C h e c k i n g #
        # ~#~#~#~#~#~#~#~ #

        self.check_if_init(["l_module_lower"],
                           sys._getframe().f_code.co_name)

        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
        # C r e a t e _ t e m pl a t e #
        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #

        l_template = []
        for m in self.l_module_lower:
            m_cap = m.capitalize()
            l_template += ["         | {0} ->".format(m_cap),
                           "           f {0}.(read, to_rst)".format(m_cap)]

        # ~#~#~#~#~#~ #
        # R e t u r n #
        # ~#~#~#~#~#~ #

        return "\n".join(l_template)

    def create_qp_write(self):
        """
        Generate write
        """
        # ~#~#~#~#~#~#~#~ #
        # C h e c k i n g #
        # ~#~#~#~#~#~#~#~ #

        self.check_if_init(["l_module_lower"],
                           sys._getframe().f_code.co_name)

        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
        # C r e a t e _ t e m pl a t e #
        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #

        str_ = "        | {0} -> write {0}.(of_rst, write) s"
        l_template = [str_.format(m.capitalize()) for m in self.l_module_lower]

        # ~#~#~#~#~#~ #
        # R e t u r n #
        # ~#~#~#~#~#~ #

        return "\n".join(l_template)

    def create_qp_tasks(self):
        """
        Generate taks
        """
        # ~#~#~#~#~#~#~#~ #
        # C h e c k i n g #
        # ~#~#~#~#~#~#~#~ #

        self.check_if_init(["l_module_lower"],
                           sys._getframe().f_code.co_name)

        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #
        # C r e a t e _ t e m pl a t e #
        # ~#~#~#~#~#~#~#~#~#~#~#~#~#~# #

        l_template = ["      {0} ; ".format(m.capitalize())
                      for m in self.l_module_lower]

        # ~#~#~#~#~#~ #
        # R e t u r n #
        # ~#~#~#~#~#~ #

        return "\n".join(l_template)
