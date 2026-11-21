====
json
====

JSON output of the results of the calculation.

The provider :c:data:`json_filename` is the next free numbered file of
the `json` directory of the |EZFIO| database, and
:c:data:`json_unit` is the unit number of that file, created open for
writing. The calculation results are appended there with the format
providers: :c:data:`json_int_fmt`, :c:data:`json_real_fmt`,
:c:data:`json_str_fmt`, :c:data:`json_true_fmt`,
:c:data:`json_false_fmt`, and the open/close variants for arrays and
dictionaries. The variants suffixed with an `x` end the line without
a trailing comma, to mark the last entry. The subroutine
:c:func:`json_close` writes the closing brace and closes the file.
