static char sccsid_sitohex__c[] =
    "@(#)sitohex_.c	1.2 96/04/09 19:13:32";

/* sitohex(in,h) ---------------------------------------

        This routine converts a short integer value in "in" into 
	an ascii hex character string in the character array h.
	NOTE: h must be 5 bytes long, because c will place the
	string terminator character '\0' in the 5th byte.
*/

sitohex_(in,h)
	short int *in;
	char *h;
{
	sprintf(h,"%.4X",*in & 0xffff);
}
