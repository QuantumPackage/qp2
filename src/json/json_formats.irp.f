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
&BEGIN_PROVIDER [ character*(64), json_array_open_fmt ]
&BEGIN_PROVIDER [ character*(64), json_array_uopen_fmt ]
&BEGIN_PROVIDER [ character*(64), json_array_close_fmt ]
&BEGIN_PROVIDER [ character*(64), json_array_close_uopen_fmt ]
&BEGIN_PROVIDER [ character*(64), json_array_close_fmtx ]
&BEGIN_PROVIDER [ character*(64), json_dict_open_fmt ]
&BEGIN_PROVIDER [ character*(64), json_dict_uopen_fmt ]
&BEGIN_PROVIDER [ character*(64), json_dict_close_uopen_fmt ]
&BEGIN_PROVIDER [ character*(64), json_dict_close_fmt ]
&BEGIN_PROVIDER [ character*(64), json_dict_close_fmtx ]
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
 json_array_open_fmt   = '(''   "'',A,''": ['')'
 json_array_uopen_fmt  = '(''   ['')'
 json_array_close_fmt  = '(''   ],'')'
 json_array_close_uopen_fmt  = '(''   ], ['')'
 json_array_close_fmtx = '(''   ]'')'
 json_dict_open_fmt    = '(''   "'',A,''": {'')'
 json_dict_uopen_fmt   = '(''   {'')'
 json_dict_close_fmt   = '(''   },'')'
 json_dict_close_uopen_fmt   = '(''   }, {'')'
 json_dict_close_fmtx  = '(''   }'')'
END_PROVIDER
