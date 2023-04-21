 BEGIN_PROVIDER [ character*(64), json_int_fmt ]
&BEGIN_PROVIDER [ character*(64), json_int_fmtx ]
&BEGIN_PROVIDER [ character*(64), json_real_fmt ]
&BEGIN_PROVIDER [ character*(64), json_real_fmtx ]
&BEGIN_PROVIDER [ character*(64), json_str_fmt ]
&BEGIN_PROVIDER [ character*(64), json_str_fmtx ]
&BEGIN_PROVIDER [ character*(64), json_true_fmt ]
&BEGIN_PROVIDER [ character*(64), json_true_fmtx ]
&BEGIN_PROVIDER [ character*(64), json_false_fmt ]
&BEGIN_PROVIDER [ character*(64), json_false_fmtx ]
 implicit none
 BEGIN_DOC
 ! Formats for JSON output.
 ! x: used to mark the last write (no comma)
 END_DOC
 json_int_fmt    = '(''   "'',A,''": '',I10,'','')'
 json_int_fmtx   = '(''   "'',A,''": '',I10)'
 json_real_fmt   = '(''   "'',A,''": '',E22.15,'','')'
 json_real_fmtx  = '(''   "'',A,''": '',E22.15)'
 json_str_fmt    = '(''   "'',A,''": "'',A,''",'')'
 json_str_fmtx   = '(''   "'',A,''": "'',A,''"'')'
 json_true_fmt   = '(''   "'',A,''": true,'')'
 json_true_fmtx  = '(''   "'',A,''": true'')'
 json_false_fmt  = '(''   "'',A,''": false,'')'
 json_false_fmtx = '(''   "'',A,''": false'')'
END_PROVIDER
