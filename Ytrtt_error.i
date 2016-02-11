/*
 * TRTT : iTeraTive Reconstruction for Tomography
 *
 * Ytrtt_error.i --
 *
 * Error package for TRTT.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2011, the MiTiV Team.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can use, modify
 * and/or redistribute the software under the terms of the CeCILL-C license as
 * circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited liability.
 *
 * In this respect, the user's attention is drawn to the risks associated with
 * loading, using, modifying and/or developing or reproducing the software by
 * the user in light of its specific status of free software, that may mean that
 * it is complicated to manipulate, and that also therefore means that it is
 * reserved for developers and experienced professionals having in-depth
 * computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling
 * the security of their systems and/or data to be ensured and, more generally,
 * to use and operate it in the same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 *
 *-----------------------------------------------------------------------------
 *
 * $Id$
 * $Log$
 */

/* REQUIREMENTS ============================================================== */

/* GLOBALS =================================================================== */

/* Error codes */
TRTT_ERR_NONE                    = 1n;
TRTT_ERR_BAD_NDIMS               = 2n;
TRTT_ERR_BAD_DIM_VALUE           = 3n;
TRTT_ERR_OUT_OF_RANGE_VALUE      = 4n;
TRTT_ERR_OUT_OF_RANGE_INDEX      = 5n;
TRTT_ERR_NOT_STRING              = 6n;
TRTT_ERR_NOT_TOMOBJ              = 7n;
TRTT_ERR_NOT_TOMDATA             = 8n;
TRTT_ERR_NOT_VOXELS              = 9n;
TRTT_ERR_NOT_PIXELS              = 10n;
TRTT_ERR_PROJECTOR_FAILURE       = 11n;
TRTT_ERR_UNSUPPORTED_JOB         = 12n;
TRTT_ERR_BAD_GEOMETRY            = 13n;
TRTT_ERR_NOT_HASH                = 14n;
TRTT_ERR_MISSING_PARAMETER       = 15n;
TRTT_ERR_BAD_PROJECTOR           = 16n;

extern trtt_error_info;
/* DOCUMENT trtt_error_info

   Ytrtt_error.i is module of TRTT code for managing errors. During the
   execution of the code, a STATUS hash table can be passed in argument of every
   "trtt_" function. If not, it is locally created in the function and destroyed
   when the function returns. STATUS namely stores information about the status
   of execution.

   STATUS
   |
   |__ code
   |
   |__ msg

   CODE gives a code for the status. Different error codes exist:
   TRTT_ERR_NONE                    = 1n;
   TRTT_ERR_BAD_NDIMS               = 2n;
   TRTT_ERR_BAD_DIM_VALUE           = 3n;
   TRTT_ERR_OUT_OF_RANGE_VALUE      = 4n;
   TRTT_ERR_OUT_OF_RANGE_INDEX      = 5n;
   TRTT_ERR_NOT_STRING              = 6n;
   TRTT_ERR_NOT_TOMOBJ              = 7n;
   TRTT_ERR_NOT_TOMDATA             = 8n;
   TRTT_ERR_NOT_VOXELS              = 9n;
   TRTT_ERR_NOT_PIXELS              = 10n;
   TRTT_ERR_PROJECTOR_FAILURE       = 11n;
   TRTT_ERR_UNSUPPORTED_JOB         = 12n;
   TRTT_ERR_BAD_GEOMETRY            = 13n;
   TRTT_ERR_NOT_HASH                = 14n;
   TRTT_ERR_MISSING_PARAMETER       = 15n;
   TRTT_ERR_BAD_PROJECTOR           = 16n;

   Each code is associated with a standard message which is displayed when necessary.

   TRTT_ERR_NONE                 => "No error."
   TRTT_ERR_BAD_NDIMS            => "Bad number of dimensions."
   TRTT_ERR_BAD_DIM_VALUE        => "Bad dimension value."
   TRTT_ERR_OUT_OF_RANGE_VALUE   => "Out of range value."
   TRTT_ERR_OUT_OF_RANGE_INDEX   => "Out of range index."
   TRTT_ERR_NOT_STRING           => "Not a string."
   TRTT_ERR_NOT_TOMOBJ           => "Not a tomographic object."
   TRTT_ERR_NOT_TOMDATA          => "Not a tomographic data."
   TRTT_ERR_NOT_VOXELS           => "Not voxels."
   TRTT_ERR_NOT_PIXELS           => "Not pixels."
   TRTT_ERR_PROJECTOR_FAILURE    => "Projector failure."
   TRTT_ERR_UNSUPPORTED_JOB      => "Unsupported job."
   TRTT_ERR_BAD_GEOMETRY         => "Bad geometry."
   TRTT_ERR_NOT_HASH             => "Not a hash tab."
   TRTT_ERR_MISSING_PARAMETER    => "Missing parameter."
   TRTT_ERR_BAD_PROJECTOR        => "Bad projector."

   The function "trtt_error_init" initializes a STATUS hashstab with
   TRTT_ERR_NONE. If an error occurs, use the function "trtt_error_set" to
   update STATUS with the right error CODE. The function
   "trtt_error_display_msg" is used to print the message error, according to
   CODE. If you want to precise the error message, you can set an extra
   explanation with the argument MSG in "trtt_error_set". The error message is
   automatically displayed by "trtt_error_set" if the keyword DISP is 1. To get
   the STATUS, use the function "trtt_error_query" which returns 1 if an error
   previously occured and STATUS was consequently updated.

   SEE ALSO: trtt_error_init, trtt_error_set, trtt_error_query,
             trtt_error_display_msg.   
*/

TRTT_ERR_MESSAGES = ["No error.",
                     "Bad number of dimensions.",
                     "Bad dimension value.",
                     "Out of range value.",
                     "Out of range index.",
                     "Not a string.",
                     "Not a tomographic object.",
                     "Not a tomographic data.",
                     "Not voxels.",
                     "Not pixels.",
                     "Projector failure.",
                     "Unsupported job.",
                     "Bad geometry.",
                     "Not a hash tab.",
                     "Missing parameter.",
                     "Bad projector."];

/* USER FUNCTIONS ============================================================ */

func trtt_error_init(void)
/* DOCUMENT trtt_error_init(void)
   
   Initialize a STATUS hashtab to report errors. It is initially filled with the
   code TRTT_ERR_NONE, and no extra MSG. STATUS must be passed in argument of
   any "trtt_" function to control its execution. In case of error, it is set to
   error status CODE, which is associated with a standard message. An extra MSG
   can be given (see trtt_error_set).

   SEE ALSO: trtt_error_info, trtt_error_set.
 */
{
    status = h_new();
    trtt_error_set, status, TRTT_ERR_NONE;
    return status;
}

func trtt_error_clear(status)
/* DOCUMENT trtt_error_clear, status
   
   Clear STATUS hashtab. It is re-initialized to TRTT_ERR_NONE status.

   SEE ALSO: trtt_error_init.
 */
{
    h_set, status, code = TRTT_ERR_NONE;
    h_set, status, msg = void;
}   

func trtt_error_set(status, code, msg=, disp=)
/* DOCUMENT trtt_error_set, status, code, msg=, disp=

   Update STATUS hashtab when an error occurs in the execution of a "trtt_"
   function. CODE specifies the type of error. MSG is an extra keyword to give
   precisions in the error messages to be displayed. If DISP=1, the error
   message is directly printed.

   SEE ALSO: trtt_error_info, trtt_error_display_msg.
 */
{      
    h_set, status, code = code;
    
    if (!is_void(msg)) {
        h_set, status, msg = msg;
    }

    if (disp) {
        trtt_error_display_msg, status=status;
    }
}

func trtt_error_query(status)
/* DOCUMENT trtt_error_query, status

   Check STATUS content. Return 1 if CODE is not TRTT_ERR_NONE => an error
   occured.

   SEE ALSO: trtt_error_info.
 */
{
    code = h_get(status,"code");
    if (code != TRTT_ERR_NONE) return TRTT_TRUE;
}

func trtt_error_display_msg(status=, code=, msg=, maxlen=)
/* DOCUMENT trtt_error_display_msg, status=, code=, msg=, maxlen=

   Display an error message, according to the content of STATUS. If STATUS is
   void, the message is given by CODE, eventually precised with MSG. MAXLEN
   specifies the maximum number of characters (default is 70) printed on a line
   (see _trtt_error_msg_to_user).

   SEE ALSO: trtt_error_info, trtt_error_set, _trtt_error_msg_to_user.
*/
{
    if (is_void(maxlen)) maxlen = 70;
    
    if (!is_void(status)) {
        code = h_get(status,"code");
        error_msg = code_to_msg(code);
        add_msg = h_get(status,"msg");      
    }
    else if (!is_void(code)) {
        error_msg = code_to_msg(code);
        add_msg = msg;       
    }
    write, format="---- %s ----\n", "TRTT ERROR";
    msg_to_user, error_msg, maxlen=maxlen;
    
    if (!is_void(add_msg)) {
        msg_to_user, add_msg, maxlen=maxlen;
    }
    write, "\n";
}

func trtt_error_code_to_msg(code)
/* DOCUMENT trtt_error_code_to_msg, code

   Return the standard error message associated with CODE.

   SEE ALSO: trtt_error_info, trtt_error_display_msg.
*/
{   
    return TRTT_ERR_MESSAGES(code);
}
code_to_msg = trtt_error_code_to_msg;


/* INTERNAL LOW LEVEL FUNCTIONS ============================================= */

func _trtt_error_msg_to_user(msg, maxlen=)
/* DOCUMENT _trtt_msg_to_user, msg, maxlen=
   
   Display a message MSG to the user. It writes the string MSG, skipping a line
   when the number of characters exceeds MAXLEN. Default value for MAXLEN is 70.
*/
{
    if (is_void(maxlen)) maxlen = 70;
       
    if (strlen(msg) < maxlen) {
        write, format="%s\n", msg;
        return;
    }
        
    nmsg = "";
    tmsg = msg;

    while (strlen(nmsg) < maxlen){
        msgtok= strtok(tmsg);
        nmsg += msgtok(1)+" ";
        tmsg = msgtok(2);
    }

    write, format="%s\n", nmsg;

    if(tmsg) msg_to_user, tmsg, maxlen=maxlen;
}
msg_to_user = _trtt_error_msg_to_user;
