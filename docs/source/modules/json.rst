.. _module_json: 
 
.. program:: json 
 
.. default-role:: option 
 
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
 
 
 
EZFIO parameters 
---------------- 
 
.. option:: empty
 
    Needed to create the json directory
 
 
 
Providers 
--------- 
 
.. c:var:: json_array_close_fmt


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_array_close_fmtx


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_array_close_uopen_fmt


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_array_open_fmt


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_array_uopen_fmt


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_dict_close_fmt


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_dict_close_fmtx


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_dict_close_uopen_fmt


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_dict_open_fmt


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_dict_uopen_fmt


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_false_fmt


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_false_fmtx


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_filename


    File : :file:`json/json.irp.f`

    .. code:: fortran

        character*(128)	:: json_filename	


    Fortran unit of the JSON file

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`file_lock`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`json_unit`

 
.. c:var:: json_int_fmt


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_int_fmtx


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_real_fmt


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_real_fmtx


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_str_fmt


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_str_fmtx


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_true_fmt


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_true_fmtx


    File : :file:`json/json_formats.irp.f`

    .. code:: fortran

        character*(64)	:: json_int_fmt	
        character*(64)	:: json_int_fmtx	
        character*(64)	:: json_real_fmt	
        character*(64)	:: json_real_fmtx	
        character*(64)	:: json_str_fmt	
        character*(64)	:: json_str_fmtx	
        character*(64)	:: json_true_fmt	
        character*(64)	:: json_true_fmtx	
        character*(64)	:: json_false_fmt	
        character*(64)	:: json_false_fmtx	
        character*(64)	:: json_array_open_fmt	
        character*(64)	:: json_array_uopen_fmt	
        character*(64)	:: json_array_close_fmt	
        character*(64)	:: json_array_close_uopen_fmt	
        character*(64)	:: json_array_close_fmtx	
        character*(64)	:: json_dict_open_fmt	
        character*(64)	:: json_dict_uopen_fmt	
        character*(64)	:: json_dict_close_uopen_fmt	
        character*(64)	:: json_dict_close_fmt	
        character*(64)	:: json_dict_close_fmtx	


    Formats for JSON output.
    x: used to mark the last write (no comma)


 
.. c:var:: json_unit


    File : :file:`json/json.irp.f`

    .. code:: fortran

        integer	:: json_unit	


    Unit file for JSON output

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`file_lock`
       * :c:data:`json_filename`


 
 
Subroutines / functions 
----------------------- 
 
.. c:macro:: json_close:


    File : :file:`json/json.irp.f`

    .. code:: fortran

        subroutine json_close



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`json_unit`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`lock_io`
       * :c:func:`unlock_io`

