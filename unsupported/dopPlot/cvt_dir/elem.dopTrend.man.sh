###############################################
# Element script                              #
#                                             #
#   This script was generated by a ClearCase  #
#  clearcvt_<cvtype> program.  It is intended #
#  to be called from the cvt_script in the    #
#  same directory.                            #
#  It contains all the info and commands to   #
#  re-create this element in a ClearCase VOB  #
#                                             #
###############################################

if [ -n "$CVT_DEBUG" ] ; then set -x;fi

PPID=${1};echo $$>$CVT_TEMP_DIRECTORY/cvtpid${PPID}

abort() {
#cleanup temp files

echo aborting... 1>&2 

rm -f $CVT_TEMP_DIRECTORY/ccaseimp$$* $CVT_TEMP_DIRECTORY/cvtpid${PPID}

$SEND_COMMAND sh exit 1
exit 1
}

trap 'abort' 1 2 3 15
if [ ! -r $CVT_VOB_DIRECTORY/'dopTrend.man'@@ ]; then  # copy test

if echo cp $CVT_FROM_DIRECTORY/'dopTrend.man' SPECIAL_SUBSTITUTION|grep SPECIAL_SUBSTITUTION >/dev/null; then

  EXTRACTCMD=`echo cp $CVT_FROM_DIRECTORY/'dopTrend.man' SPECIAL_SUBSTITUTION| sed s?SPECIAL_SUBSTITUTION?$CVT_TEMP_DIRECTORY/ccaseimp$$?g`

$EXTRACTCMD

ret=$?
if [ "$ret" != 0 ]
then
	 abort
fi

else

cp $CVT_FROM_DIRECTORY/'dopTrend.man' SPECIAL_SUBSTITUTION > $CVT_TEMP_DIRECTORY/ccaseimp$$

ret=$?
if [ "$ret" != 0 ]
then
	 abort
fi

fi

 cp $CVT_TEMP_DIRECTORY/ccaseimp$$ $CVT_VOB_DIRECTORY/'dopTrend.man' 


ret=$?
if [ "$ret" != 0 ]
then
	 abort
fi

fi # copy test

if [ -d $CVT_VOB_DIRECTORY/'dopTrend.man'@@ ] ; then  # elem existence check

  output=`cleartool lshistory -since 16-Jun-98.12:23:44 -branch ${CVT_USING_CSPEC-main} $CVT_VOB_DIRECTORY/'dopTrend.man'|grep "create version"|grep -v "/0"|grep -v Destroy`

  if [ -n "$CVT_UPDATE" ] ; then
    output=
  fi

  if [ -n "$output" ] ; then

    echo Already converted 'dopTrend.man' 1>&2 
    rm -f $CVT_TEMP_DIRECTORY/ccaseimp$$* ; exit 0

  else
    INCREMENTAL=TRUE
export INCREMENTAL
  fi

    echo Already created $CVT_VOB_DIRECTORY/'dopTrend.man' 1>&2 

    cleartool uncheckout -rm $CVT_VOB_DIRECTORY/'dopTrend.man' >/dev/null 2>/dev/null

else  # elem existence check

if [ ! -r $CVT_VOB_DIRECTORY/'dopTrend.man' ] ; then touch $CVT_VOB_DIRECTORY/'dopTrend.man' ;fi 

$SEND_COMMAND "setevent -user 'olawlor' -group 'ucsl' -name 'Orion Lawlor' -host 'sparc1k' -time '16-Jun-98.12:23:44'"
$SEND_COMMAND mkelem -nc -nco -rm $CVT_VOB_DIRECTORY/'dopTrend.man'
$SEND_COMMAND protect -chmod 775 $CVT_VOB_DIRECTORY/'dopTrend.man'@@

$SEND_COMMAND sh cp /dev/null ${MK_VER_FILE}

while [ ! -r "${MK_VER_FILE}" ]
 do
   if echo; then :;else echo 1>&2 elem script exiting;abort;fi
   sleep 1
 done


rm -f ${MK_VER_FILE}
:

fi # elem exist check

rm -f $MK_VER_FILE


$SEND_COMMAND "setevent -user 'olawlor' -group 'ucsl' -name 'Orion Lawlor' -host 'sparc1k' -time '16-Jun-98.12:23:44'"

EXTRACTCMD=`echo cp $CVT_FROM_DIRECTORY/'dopTrend.man' SPECIAL_SUBSTITUTION`
export EXTRACTCMD

cvt_mkver.sh $CVT_VOB_DIRECTORY/'dopTrend.man' 16-Jun-98.12:23:44 main "$CVT_TEMP_DIRECTORY/ccaseimp$$" "one"  <<"_COMMENT_"
made from unix file
.
_COMMENT_

ret=$?
if [ "$ret" != 0 ]
then
	 abort
fi


while [ ! -r "${MK_VER_FILE}" ]
 do
   if echo; then :;else echo 1>&2 elem script exiting;abort;fi
   sleep 1
 done


if rm $MK_VER_FILE ; then :
else echo 1>&2 Unexpected error;abort;fi

rm -f $CVT_TEMP_DIRECTORY/ccaseimp$$*


