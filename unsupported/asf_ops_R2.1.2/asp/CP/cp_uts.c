/* Alaska SAR Processor (ASP) %W% %E% %U% */

/* CV 3/6/96 - Add the flag "ucs_product" to process the UNCOMPENSATED */
/*		product */
/* CV 2/96 - fill in MODE, DATASET, VERSION in PMF file */
/*	   - Updated the format of scan_result_file */ 
/*	   - Correct the nrec & nreclen of image file in PMF file */
/*	     (add 1 into nrec & 192 bytes into nreclen */
/* CV 1/96 - Modify the CLEANUP_REQUEST */ 	
/*	   - Updated the ODL Object names with new DD v1.3 */
/* CLW 12/18/95 - Patch more PMF paramters from system engineer */
/* CLW 12/14/95 - Extract version from ceos_ldr filename */
/* CLW 12/14/95 - Add CAL_POW extraction from scan result file */
/* CLW 12/13/95 - Add CLEANUP_REQUEST function */
/* CLW 12/08/95 - Replace ASP_REQUEST_ACK & ASP_REQUEST_STATUS to
		  SUBSYSTEM_ACK & SUBSYSTEM_STATUS requested by CP */
/* CLW 12/06/95 - Modify product names due to changes made by SE */
/* CLW 12/05/95 - Modify to accept radarsat mode WD1, WD2, WD3 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <syslog.h>
#include <math.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include "asf.h"
#include "cp_uts.h"
#include <pthread_wrapper.h>
#include <procdec.h>

#define SP 32
#define FAIL -1
#define PASS 0

extern int vbose;
extern int chkrds;
extern char ASP_VERSION[8];
extern RQST_PTR Cur_Rqst;              /* current request */
extern SV sv1,sv2;                     /* state vectors */
extern TC_FILE tcf;                    /* time correlation element */
extern char PROC_PATH[];
extern char JOBS_PATH[];
extern int ucs_product;
extern char chk_media_id[];
extern int Write_CP_Error_Status;
extern int st_fmt;
extern int asp_job_id, asp_rev_id;
extern char job_name[100];
extern char scan_result_version[30];
extern char scan_result_filename[50];
extern char cal_params_filename[50];
extern char cal_status[50];
extern char cal_comment[300];
extern double noise_fct;
extern double linear_conv_fct;
extern double offset_conv_fct;
extern double islr_rng, pslr_rng;
extern double range_res, azimuth_res;
extern double range_ambig, azimuth_ambig;
extern double ori_error, dist_skew;
extern double crosst_scale, alongt_scale;
extern double crosst_locerr, alongt_locerr;
extern double rel_radio_acc;
extern char frame_mode[];  /* amm */
extern char host[16];

/* ASP-CP 2.1, 11/21/96 */
extern char data_direction[];
extern char media_type[];
extern int  experimental_beam;
extern char compensation_flag[];
extern char quicklook_flag[];
extern char inlog[]; /* defined in tape_mnt.c */
extern int use_rds_scan;
extern int rds_center_format;
extern int max_agc;

static CP_RQST_PTR cp_rp;
static CAL_POW cal;
static char *save_rqst;
static char *save_seg;
static char *save_frame, *src;
static int asp_cp = 1;
static int quicklook = 0;
static char *asp_stage = "IDLE";
static char rev_name[20];
int first_st_blk, first_st_bit;

/* handle_proc_req() -------------------------------------------------------
        This routine is triggered by cp job request signal.
        It will performs:
                - generate job file on disk
                - process job
*/
void handle_proc_req( ODL msg )
{
	ODL ack;
	int status;
	char job_fn[132];
	FILE *fp;

	ucs_product = 0;
	quicklook = 0;
	save_rqst = (char *) ODLToStr(msg,NULL);
	if ( asp_cp ){
		asp_msg(LOG_DEBUG,"Acknowledge CP...");
		if ((ack = GetAckMsg(msg, "SUBSYSTEM_ACK")) == NULL) {
			asp_exit(CP_ACK_ERR);
		}
		WriteMsgToServer(ack);
		ODLFree(ack);
	}
	strcpy(asp_stage,"IDLE");
	asp_msg(LOG_DEBUG,"Handling cp processing request...");

	if ( cp_rp->proc_type == FRAME_RQST ) 
		status = GetFrameRqst(msg,cp_rp,Cur_Rqst,&sv1,&tcf);
	else if ( cp_rp->proc_type == SCAN_RQST )
		status = GetScanRqst(msg,cp_rp,Cur_Rqst,&sv1,&tcf);

	if ( status == FAIL ){
/* ASP-CP 2.1, 11/21/96 */
             if (Write_CP_Error_Status == FAIL)   
		asp_msg(LOG_ERR,asp_messages[CP_RQST_ERR]);
             else
                 Write_CP_Error_Status = FAIL;
		if ( asp_cp ) 
		   if (WriteStatusToCP(msg, status)==FAIL)
			asp_msg(LOG_ERR,asp_messages[WRITE_CPSTAT_ERR]);
		strcpy(asp_stage,"IDLE");
		return;
	} else {
		sprintf(job_fn,"%s%s.odl",JOBS_PATH,Cur_Rqst->jobname);
        	if ( (fp = fopen(job_fn,"w+"))==NULL)
	   		printf("Can't open %s file\n",job_fn);
        	else fprintf(fp,"%s",save_rqst);
		fclose(fp);
		write_job_file( Cur_Rqst, &sv1, &tcf );	
	}
	if ( cp_rp->proc_type )
		strcpy(asp_stage,"PROCESS");
	else
		strcpy(asp_stage,"SCAN");

        status = process_cp_job();
        if (status == FAIL)  status = Write_CP_Error_Status;
	if ( asp_cp ) { 
		if (WriteStatusToCP(msg, status)==FAIL)
			asp_msg(LOG_ERR,asp_messages[WRITE_CPSTAT_ERR]);
        }
     
        Write_CP_Error_Status = FAIL;
	strcpy(asp_stage,"IDLE");

} /* end of handle_proc_req */

void handle_frame_req( ODL msg, void *arg )
{

	asp_signals();
	asp_msg(LOG_DEBUG,"ASP received a frame request.");
    	cp_rp = (CP_RQST_PTR) malloc(sizeof(CP_RQST_RCD));
	Cur_Rqst = &cp_rp->rp;
	cp_rp->proc_type = FRAME_RQST;
	handle_proc_req(msg);
	free(cp_rp);

} /* end of handle_frame_req */

void handle_scan_req( ODL msg, void *arg )
{
	float x;

	asp_signals();
	asp_msg(LOG_DEBUG,"ASP received a scan request.");
    	cp_rp = (CP_RQST_PTR) malloc(sizeof(CP_RQST_RCD));
	Cur_Rqst = &cp_rp->rp;
	cp_rp->proc_type = SCAN_RQST;
	handle_proc_req(msg);
	free(cp_rp);

} /* end of handle_scan_req */

void handle_heartbeat( ODL msg, void *arg )
{
	ODL ack;
	int i;

        for ( i = 0; (asp_proc_stage[i] != NULL) &&
                        strcmp(asp_stage,asp_proc_stage[i]); ++i )
	  ;

	printf("asp_proc=%s\n",asp_proc_stage[i]);
        if ((ack = GetAckMsg(msg, "SUBSYSTEM_STATUS")) == NULL) {
		asp_msg(LOG_DEBUG,"can't get SUBSYSTEM_STATUS\n");
        }
	if ( !ODLSetVal(ack, "SUBSYSTEM_STATUS.BODY.SUB_TYPE", "HEARTBEAT") ){
		asp_msg(LOG_DEBUG,"can't set BODY.SUB_TYPE\n");
	}
/* ASP-CP 2.1, 11/21/96 */
	if ( !ODLSetVal(ack, "SUBSYSTEM_STATUS.BODY.STATUS", "COMPLETED") ){
/*
	if ( !ODLSetVal(ack, "SUBSYSTEM_STATUS.BODY.STATUS", PASS) ){
*/
		asp_msg(LOG_DEBUG,"can't set BODY.STATUS\n");
	}
	if ( !ODLSetVal(ack, "SUBSYSTEM_STATUS.BODY.COMMENT",asp_proc_stage[i]) ){
		asp_msg(LOG_DEBUG,"can't set BODY.COMMENT\n");
	} 
	WriteMsgToServer(ack);
	ODLFree(ack);

} /* end of handle_heartbeat */

void handle_cleanup_req( ODL msg, void *arg )
{
	int status = PASS;
	char *save_flag, *save_dir;
	int job_id, len, itmp, rev;
	struct dirent *dp;
	DIR *dirp, *opendir();
	char *jobname = NULL;
	char job_rev[20], *string, rev_temp[10];

	asp_signals();
	asp_msg(LOG_DEBUG,"ASP received a cleanup request.");
	msg = Value(msg);
        string = ODLGetString(msg, SATELLITE, &status);
	if (vbose) printf("SATELLITE= %s, ", string);
        if ( status == FAIL ) 
          asp_msg(LOG_DEBUG,"can't get %s", SATELLITE);
        else {
          sprintf(job_rev,"%.1s",&string[0]);
          rev = ODLGetInt(msg, REV, &status);
	  if (vbose) printf("rev=%d, ", rev);
	  sprintf (rev_temp,"%0.5d",rev);
          if ( status == FAIL ) asp_msg(LOG_DEBUG,"can't get %s", REV);
	  else {
            strcat(job_rev,rev_temp);
            if (vbose) printf("job_rev=%s\n", job_rev);
            job_id = ODLGetInt(msg, PROC_REQ_ID, &status);
            if (status == FAIL) 
		asp_msg(LOG_DEBUG,"can't get %s", PROC_REQ_ID);
            else {
              dirp = opendir(JOBS_PATH);
              while ((dp = readdir(dirp)) != NULL){
                len = strlen(dp->d_name);
                if ( len != 14 ) continue;
                if (strncmp(dp->d_name, job_rev,6)) continue;
                else {
                  itmp = atoi(&dp->d_name[9]);
                  if ( itmp != job_id ) continue;
                  jobname = (char *) malloc(15);
                  strcpy(jobname,dp->d_name);
                  if ( vbose ) printf("Found jobname %s\n",jobname);
                }
              }
	      closedir(dirp);
              if ( jobname == NULL){
		if (!strcmp(asp_stage,"IDLE")) status = PASS; 
                else {
		    asp_msg(LOG_DEBUG,"Unable to find REV %s, JOB_ID %d",
                        job_rev, job_id);
                    status = FAIL;
		}
              } else {
                save_flag = ODLGetString(msg, SAVE_FLAG, &status);
                if ( status == FAIL ){
                   asp_msg(LOG_DEBUG,"can't get %s", SAVE_FLAG);
		   free(save_flag);
		}
                else {
                   if ( !strcmp(save_flag,"NO") ){
			free(save_flag);
                        cleanup_data(jobname);
                   } else {
		       free(save_flag);
                       save_dir = ODLGetString(msg, SAVE_DIR, &status);
                       if ( status == FAIL ) {
                           asp_msg(LOG_DEBUG,"can't get %s", SAVE_DIR);
			   free(save_dir);
			}
                        else {
			   dirp = opendir(save_dir);
        		   if ((dp = readdir(dirp)) == NULL){
           			asp_msg(LOG_DEBUG,"can not find the SAVE_DIR:%s"					, save_dir);
				free(save_dir);
				status = FAIL;
			   } else
                                save_data( jobname, save_dir );
			  free(save_dir);
			  closedir(dirp);
                        }
                   }
                }
              }
	    }
          }
       }
       if ( status == FAIL )
		asp_msg(LOG_ERR,asp_messages[INVALID_CLEANUP_RQST]);
       if ( asp_cp ) {
		if (WriteCleanupStatusToCP(msg, status, job_id)==FAIL)
			asp_msg(LOG_ERR,asp_messages[WRITE_CPSTAT_ERR]);
		else asp_msg(LOG_DEBUG,"ASP finish cleanup/save the job ");
       }

} /* end of handle_cleanup_req */

int activate_cp(cp_argc,cp_argv)
	int cp_argc;
	char *cp_argv[];
{
	AsfApp app;
	register int i;
	ODL msg, GetODLFile();
	char *RqstType;
	char RqstFileName[80];
	int status;

	for (i=1; cp_argv[i] != NULL; ++i){
		if ( vbose )
			printf(">%s< %d ",cp_argv[i],strlen(cp_argv[i]));
		if ( !strcmp(cp_argv[i], "-asp_cp") ) {
			asp_cp = 0;
			sprintf(RqstFileName,"%s%s",JOBS_PATH,cp_argv[i+1]);
		}
	}
	app = AsfAppInitialize(ASP_VERSION, &cp_argc, cp_argv);
        if ( app == NULL ){
	    if ( asp_cp ) {
                asp_msg(LOG_ERR,asp_messages[ASF_INIT_ERR]);
		return(FAIL);
	    }
        }
	if ( asp_init() == FAIL ) return(FAIL);
	if ( asp_cp ){
        	AsfAddCallback( app, FRAME_REQUEST, handle_frame_req,
                                NULL, 0 );
        	AsfAddCallback( app, SCAN_REQUEST, handle_scan_req,
                                NULL, 0 );
        	AsfAddCallback( app, CLEANUP_REQUEST, handle_cleanup_req,
                                NULL, 0 );
        	AsfAddCallback( app, SUBSYSTEM_HEARTBEAT, handle_heartbeat,
                                NULL, 0 );
        	AsfAppMainLoop(app);
	} else {
		msg = (ODL) GetODLFile( RqstFileName );	
    		RqstType = ODLGetString(msg,FRAME_TYPE, &status);
		if (status == FAIL) {
    		     RqstType = ODLGetString(msg,SCAN_TYPE, &status);
		     if (status == FAIL){
			free(RqstType);
			asp_msg(LOG_DEBUG,"can't get %s", "MSG_TYPE"); 
		     }
		}
		if ( !strcmp(RqstType,FRAME_REQUEST) ){
			free(RqstType);
			handle_frame_req( msg, NULL );
		}
		else if ( !strcmp(RqstType,SCAN_REQUEST) ){
			free(RqstType);
			handle_scan_req( msg, NULL );
		}
		/*ODLFree(msg);*/
		return(PASS);
	}
}
int GetODLRqst (rqst, crp, rp, sv, tcf )
 	ODL rqst;
	CP_RQST_PTR crp;
	RQST_PTR rp;
	SV_PTR sv;
	TC_FILE_PTR tcf;
{
	char *string;
	int i, itmp, month, day;
	int is, fs;
	static char *mtxt[] = { "xxx", "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
				"JUL", "AUG", "SEP", "OCT", "NOV", "DEC" };
	double angle;
	char gha_fn[132];
	FILE *fp;
	GMT gmt;
	int status;
	double deltat, get_gmt_diff();

/* ------------------------ Fill in rp data block ---------------- */

/* Get take_id; format is "SS/X/rrrrr.nn" */
    string = ODLGetString(rqst, SATELLITE, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", SATELLITE); 
	free(string);
	return(FAIL);
    } else sprintf( rp->take_id, "%s/", string );
    free(string);

    string = ODLGetString(rqst, SENSOR, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", SENSOR); 
	free(string);
	return(FAIL);
    } else {
	strcat( rp->take_id, string );
	strcat( rp->take_id, "/" );
    }

    itmp = ODLGetInt( rqst, REV, &status );
    if (status == FAIL) {
	asp_msg(LOG_DEBUG,"can't get %s", REV); 
	return(FAIL);
    } else {
	sprintf(string,"%05.5d", itmp);
	strcat( rp->take_id, string );
	strcat( rp->take_id, "." );
    }

    itmp = ODLGetInt( rqst, SEQ, &status );
    if (status == FAIL) {
	asp_msg(LOG_DEBUG,"can't get %s", SEQ); 
	free(string);
	return(FAIL);
    } else {
	sprintf(string,"%02.2d",itmp);
	strcat( rp->take_id, string );
    }
    free(string);

    if (vbose) printf("Got rp->take_id %s\n",rp->take_id);

/* Get tape_id; format is "WSnnnn" */

    string = ODLGetString(rqst, MEDIA_ID, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", MEDIA_ID); 
        free(string);
	return(FAIL);
    } else {
            strcpy( rp->tape_id, string );
/*
        if (strlen(string) <= 7)
              strcpy( rp->tape_id, string );
        else {
              strncpy( rp->tape_id, string+2, 2 ); /*WS
              strcpy( rp->tape_id+2, string+8); /*xxxx
       }
*/
    }
    free(string);
    if (vbose) printf("Got rp->tape_id %s\n",rp->tape_id);

/* ASP-CP 2.1, 11/21/96 */

    string = ODLGetString(rqst, MEDIA_TYPE, &status);  /* DCRSI/ID-1 */
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", MEDIA_TYPE); 
        free(string);
	return(FAIL);
    } else strcpy( media_type, string );
    free(string);
    if (vbose) printf("Got MEDIA_TYPE %s\n", media_type);

    string = ODLGetString(rqst, DATA_DIRECTION, &status); /*FORWARD/REVERSE*/
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", DATA_DIRECTION); 
        free(string);
	return(FAIL);
    } else strcpy( data_direction, string );
    free(string);
    if (vbose) printf("Got DATA_DIRECTION %s\n",data_direction);

/* add check_media_id gwu 9/20/96 */

    string = ODLGetString(rqst, CHECK_MEDIA_ID, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", CHECK_MEDIA_ID); 
        free(string);
	return(FAIL);
    } else strcpy( chk_media_id, string );
    if (vbose) printf("Got CHECK_MEDIA_ID %s\n",chk_media_id);

    rp->start_blk = ODLGetInt( rqst, START_BLOCK, &status );
    if (status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", START_BLOCK); 
        free(string);
	return(FAIL);
    }; 
    if (vbose) printf("Got rp->start_blk %d\n",rp->start_blk);
    rp->end_blk = ODLGetInt( rqst, END_BLOCK, &status );
    if (status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s",  END_BLOCK); 
        free(string);
	return(FAIL);
    }
    if (vbose) printf("Got rp->end_blk %d\n",rp->end_blk);

/* Get start & end time */

    if (GetODLGMT(rqst, START_TIME, &rp->start) == FAIL ) return(FAIL);

    if (GetODLGMT(rqst, END_TIME, &rp->end) == FAIL ) return(FAIL);
	
/*			  12345678901234567890123456789  */
/* Get job id; format is "Srrrrrnn xyyyyy ddmmmyyyy zzz" */
/* Only the first 15 characters are used */
/* ??? Assume yyyy is PROC_REQ_ID */
    
    crp->job_id = ODLGetInt(rqst, PROC_REQ_ID, &status);
    if ( status == FAIL ) {
	asp_msg(LOG_DEBUG,"can't get %s", PROC_REQ_ID); 
        free(string);
	return(FAIL);
    } else {
	sprintf(string,"%05.5d",crp->job_id);
    }
    get_month_day( rp->start.yr, rp->start.day, &month, &day );
    if ( vbose ) printf("month %d day %d\n", month, day);
    if ( !strcmp(rp->type,asp_product_type[CPX]) ){
    	sprintf(rp->id,"%.1s%.5s%.2s X%.5s %.2d%.3s%.4d %.3d", 
		&rp->take_id[0], &rp->take_id[5], &rp->take_id[11], 
		string, day, mtxt[month], rp->start.yr, 
		(int)rp->start.second );
    } else if ( ucs_product) {
    	sprintf(rp->id,"%.1s%.5s%.2s %.1s%.5s %.2d%.3s%.4d %.3d", 
		&rp->take_id[0], &rp->take_id[5], &rp->take_id[11], 
		"N", string, day, mtxt[month], rp->start.yr,
		(int)rp->start.second );
    } else if ( quicklook ) {
    	sprintf(rp->id,"%.1s%.5s%.2s %.1s%.5s %.2d%.3s%.4d %.3d", 
		&rp->take_id[0], &rp->take_id[5], &rp->take_id[11], 
		"Q", string, day, mtxt[month], rp->start.yr,
		(int)rp->start.second );
    } else {
    	sprintf(rp->id,"%.1s%.5s%.2s %.1s%.5s %.2d%.3s%.4d %.3d", 
		&rp->take_id[0], &rp->take_id[5], &rp->take_id[11], 
		&rp->type[0],
		string, day, mtxt[month], rp->start.yr,
		(int)rp->start.second );
    }
    free(string);
    sprintf(rp->jobname,"%8.8s%6.6s",&rp->id[0],&rp->id[9]);
    rp->jobname[8] = tolower( rp->jobname[8] );
    if ( vbose ) printf("Got rp->id %s\n", rp->id);
    if ( vbose ) printf("Got rp->jobname %s\n", rp->jobname);

    crp->mode = ODLGetString(rqst, JOB_MODE, &status);
    if (vbose) printf("mode (string) =%s\n",crp->mode);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", JOB_MODE); 
	return(FAIL);
    }

/* ASP-CP 2.1, 11/21/96 */
    experimental_beam = ODLGetInt(rqst, EXPERIMENTAL_BEAM, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", EXPERIMENTAL_BEAM); 
	return(FAIL);
    }
    if (vbose) printf("EXPERIMENTAL_BEAM = %d\n", experimental_beam);

/* ------------------------ Fill in sv data block ---------------- */

    if (GetODLGMT(rqst, SV_TIME, &sv->gmt) == FAIL ) return(FAIL);
    else {
	deltat = get_gmt_diff(&rp->start,&sv->gmt);
	if ( deltat > 6000. ){ 
		asp_msg(LOG_ERR,asp_messages[RQST_SV_ERR]);
		asp_msg(LOG_DEBUG,"sv & rqst start time difference: %f",deltat);
		return(FAIL);
	}
    }
    sv->rev = ODLGetInt(rqst, SV_REV, &status);
    if ( status == FAIL ) {
	asp_msg(LOG_DEBUG,"can't get %s", SV_REV); 
	return(FAIL);
    }
    if (vbose) printf("Got sv->rev %d\n", sv->rev);
    sv->pos.x = ODLGetDouble(rqst,SV_X_POS, &status);
    if ( status == FAIL ) {
	asp_msg(LOG_DEBUG,"can't get %s", SV_X_POS); 
	return(FAIL);
    }
    sv->pos.y = ODLGetDouble(rqst, SV_Y_POS, &status);
    if ( status == FAIL ) {
	asp_msg(LOG_DEBUG,"can't get %s", SV_Y_POS); 
	return(FAIL);
    }
    sv->pos.z = ODLGetDouble(rqst, SV_Z_POS, &status);
    if ( status == FAIL ) {
	asp_msg(LOG_DEBUG,"can't get %s", SV_Z_POS); 
	return(FAIL);
    }
    if (vbose) printf("Got sv->pos %g %g %g\n", sv->pos.x, sv->pos.y,
						sv->pos.z);
    sv->vel.x = ODLGetDouble(rqst, SV_X_VEL, &status);
    if ( status == FAIL ) {
	asp_msg(LOG_DEBUG,"can't get %s", SV_X_VEL); 
	return(FAIL);
    }
    sv->vel.y = ODLGetDouble(rqst, SV_Y_VEL, &status);
    if ( status == FAIL ) {
	asp_msg(LOG_DEBUG,"can't get %s", SV_Y_VEL); 
	return(FAIL);
    }
    sv->vel.z = ODLGetDouble(rqst, SV_Z_VEL, &status);
    if ( status == FAIL ) {
	asp_msg(LOG_DEBUG,"can't get %s", SV_Z_VEL); 
	return(FAIL);
    }
    if (vbose) printf("Got sv->vel %g %g %g\n", sv->vel.x, sv->vel.y,
						sv->vel.z);

    string = ODLGetString(rqst, SV_TYPE, &status );
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", SV_TYPE); 
        free(string);
	return(FAIL);
    } else {
	for ( i = SV_PREC_LNUM; i < SV_PREC_HNUM; i++ ){
	   if ( !strcmp( string, sv_precision[i] ) ) {
		sv->precision = i;
		break;
	   }
	}
        free(string);
	if (( i == SV_PREC_HNUM ) || (i > SV_PREC_HNUM)) return(FAIL);
    }
    
/* ------------------------ Fill in tcf data block ---------------- */
    tcf->rev = ODLGetInt(rqst,TCF_REV, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", TCF_REV); 
	return(FAIL);
    } 
    if (GetODLGMT(rqst, TCF_TIME, &tcf->gmt) == FAIL ) return(FAIL);

    tcf->bt = (unsigned int) ODLGetInt(rqst, TCF_BT, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", TCF_BT); 
	return(FAIL);
    } 
    tcf->delta = ODLGetInt(rqst, TCF_DELTA, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", TCF_DELTA); 
	return(FAIL);
    } 
    if (vbose)
	printf("tcf: rev %d bt %u delta %d\n", tcf->rev,
			tcf->bt, tcf->delta );


/* ----------- Get input/output filenames */

    crp->frame_file = ODLGetString(rqst, SCAN_RESULTS_FILE, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", SCAN_RESULTS_FILE); 
	return(FAIL);
    }
    if (vbose)
	printf("frame_file %s\n", crp->frame_file);

/* ------------------------ Write GHA file ---------------- */

    if (GetODLGMT(rqst, GHA_TIME, &gmt) == FAIL ) return(FAIL);

    angle = ODLGetDouble(rqst, GHA_ANGLE, &status);
    if (status == FAIL){
	asp_msg(LOG_DEBUG,"can't get %s", GHA_ANGLE);
	return(FAIL);
    } 
    string = (char *) malloc(80);
    sprintf(gha_fn,"%s%s",PROC_PATH,ASP_GHA_FILE);
    sprintf(string,"cp %s %s.sav\n", gha_fn, gha_fn);
    system(string);
    if((fp = fopen(gha_fn,"a+"))==NULL){
	if ( vbose ) printf("Can't open %s file\n",gha_fn);
	return (FAIL);
    }
    is = (int)(gmt.second);
    fs = (int)((gmt.second - (float)(is))*1000.0+0.5);
    fprintf(fp,"%5d:%03.3d:%02.2d:%02.2d.%02.2d.%03.3d   %10.6f\n",
		gmt.yr,gmt.day,gmt.hr,gmt.min,is,fs,angle);
    fclose(fp);
    sprintf(string,"sort -u %s > %sgha_tmp\n", gha_fn, PROC_PATH);
    system(string);
    sprintf(string,"mv %sgha_tmp %s\n", PROC_PATH, gha_fn);
    system(string);
    sprintf(string,"chmod 774 %s\n", gha_fn);
    system(string);
    free(string);

    if (vbose) printf("Updated %s file\n",gha_fn);

/* Save additional parameters */

    string = ODLGetString(rqst, ACTIVITY, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", ACTIVITY); 
    	free(string);
	return(FAIL);
    } else {
	for ( i = 0; (cp_activity[i] != NULL) && 
			strcmp(string,cp_activity[i]); ++i );
	if ( cp_activity[i] != NULL ) crp->activity = i;
	else {
		asp_msg(LOG_DEBUG,"unsupported %s", ACTIVITY);
                printf("****ACTIVITY=%s\n", string);
    		free(string);
		return(FAIL);
	}
        free(string);
    }
    return(PASS);

} /* end of GetODLRqst */

int GetScanRqst (rqst, crp, rp, sv, tcf )
 	ODL rqst;
	CP_RQST_PTR crp;
	RQST_PTR rp;
	SV_PTR sv;
	TC_FILE_PTR tcf;
{
	int status;
	char *platform, *string;
	int asp_job_id, asp_rev_id;

/* Blanks or Zeros the rest of parameters */

    rqst = Value(rqst);
    strcpy(rp->type,"STD");
    memset(rp->site,SP,sizeof(rp->site));
    rp->ave_hght = 0.;
    rp->proc_gain = 0;
    strcpy(rp->deskew,"NOT");
    strcpy(rp->gnd_slnt_rg, "SLANT");
    rp->lat = rp->lon = -1000.; 

    asp_job_id = ODLGetInt(rqst, PROC_REQ_ID, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", PROC_REQ_ID);
	return(FAIL);
    }
    asp_rev_id = ODLGetInt(rqst, REV, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", REV);
	return(FAIL);
    }
    platform = ODLGetString(rqst, SATELLITE, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", SATELLITE);
	free(platform);
	return(FAIL);
    }
    sprintf(job_name, "REV %s %.05dS%.05d", platform,asp_rev_id,asp_job_id);
    asp_msg(LOG_DEBUG,"%s", job_name);
    free(platform);
    string = ODLGetString(rqst, FRAME_MODE, &status);
    if ( status == FAIL ){
        asp_msg(LOG_DEBUG,"can't get %s", FRAME_MODE);
        free(string);
        return(FAIL);
    } else {
        strcpy(frame_mode,string);
        free(string);
    }
    if (vbose) printf("frame_mode = %s\n",frame_mode);
    return(GetODLRqst(rqst,crp,rp,sv,tcf));

} /* end of GetScanRqst */

int GetFrameRqst (rqst, crp, rp, sv, tcf )
 	ODL rqst;
	CP_RQST_PTR crp;
	RQST_PTR rp;
	SV_PTR sv;
	TC_FILE_PTR tcf;
{
	char *string, *string_product_type, stmp[132], *platform;
	int i;
	int status;
	int asp_job_id, asp_rev_id;

    rqst = Value(rqst);
    
    asp_job_id = ODLGetInt(rqst, PROC_REQ_ID, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", PROC_REQ_ID);
	return(FAIL);
    }
    asp_rev_id = ODLGetInt(rqst, REV, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", REV);
	return(FAIL);
    }
    platform = ODLGetString(rqst, SATELLITE, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", SATELLITE);
        free(platform);
	return(FAIL);
    }
    sprintf(rev_name,"%s%.05d", platform,asp_rev_id);
    if (vbose) printf("rev_name=%s\n",rev_name);

/* amm */
    string = ODLGetString(rqst, FRAME_MODE, &status);
    if ( status == FAIL ){
        asp_msg(LOG_DEBUG,"can't get %s", FRAME_MODE);
        free(string);
        return(FAIL);
    } else {
        strcpy(frame_mode,string);
        free(string);
    }
    if (vbose) printf("frame_mode = %s\n",frame_mode);

    if (strcmp(frame_mode,"ARCTIC")){  /* add CV 5/16/97 */
	sprintf(asp_error_msg, "It is not ARCTIC frame mode");
	asp_msg(LOG_ERR,asp_error_msg);
        free(string);
        return(FAIL);
    }

/* ASP-CP 2.1, 11/21/96 */
/* changed string to "string_product_type" */

    string_product_type = ODLGetString(rqst, PRODUCT_TYPE, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", PRODUCT_TYPE); 
	free(string_product_type);
	return(FAIL);
    } else {
	for ( i = 0; (cp_product_type[i] != NULL) && 
			strcmp(string_product_type,cp_product_type[i]); ++i );
	if ( cp_product_type[i] != NULL ) 
			strcpy( rp->type, asp_product_type[i] );
	else {
		asp_msg(LOG_DEBUG,"invalid %s", PRODUCT_TYPE);
		free(string_product_type);
		return(FAIL);
	}
	free(string_product_type);
    }
    if (vbose) printf("Got rp->type %s\n",rp->type);

/* ASP-CP 2.1, 11/21/96 */

    string = ODLGetString(rqst, COMPENSATION_FLAG, &status);
    if ( status == FAIL ){
        asp_msg(LOG_DEBUG,"can't get %s", COMPENSATION_FLAG);
        free(string);
        return(FAIL);
    } else {
        strcpy(compensation_flag,string);
        free(string);
    }
    if (vbose) printf("COMPENSATION_FLAG = %s\n", compensation_flag);

    string = ODLGetString(rqst, QUICKLOOK_FLAG, &status);
    if ( status == FAIL ){
        asp_msg(LOG_DEBUG,"can't get %s", QUICKLOOK_FLAG);
        free(string);
        return(FAIL);
    } else {
        strcpy(quicklook_flag,string);
        free(string);
    }
    if (vbose) printf("QUICKLOOK_FLAG = %s\n", quicklook_flag);


    if (!strcmp(string_product_type, "RAMP")) {
       if (!strcmp(compensation_flag, "YES"))
        Write_CP_Error_Status = 25; /*  25 != FAIL */
        sprintf(asp_error_msg, "RAMP/AMM products are not supported by this version of ASP");
        asp_msg(LOG_ERR, asp_error_msg);
        free(string_product_type);
        return(FAIL);
    }

    if (!strcmp(string_product_type, "STANDARD")) {
       if (!strcmp(compensation_flag, "NO"))
           strcpy(rp->type, "NUC");
    }

    if ((!strcmp(string_product_type, "STANDARD")) |
        (!strcmp(string_product_type, "RAMP"))) {
       if (!strcmp(quicklook_flag, "YES"))
           strcpy(rp->type, "QLK");
    }

    if (!strcmp(quicklook_flag, "YES")) {
      if ((!strcmp(string_product_type, "COMPLEX")) |
          (!strcmp(string_product_type, "CCSD"))) {
        Write_CP_Error_Status = 25; /*  25 != FAIL */
        sprintf(asp_error_msg, "No QUICKLOOK for product type %s", string_product_type);
        asp_msg(LOG_ERR, asp_error_msg);
        free(string_product_type);
        return(FAIL);
      }
    }

    if (!strcmp(compensation_flag, "NO")) {
      if ((!strcmp(string_product_type, "COMPLEX")) |
          (!strcmp(string_product_type, "CCSD"))) {
        Write_CP_Error_Status = 25; /*  25 != FAIL */
        sprintf(asp_error_msg, "No UNCOMPENSATED for product type %s", string_product_type);
        asp_msg(LOG_ERR, asp_error_msg);
        free(string_product_type);
        return(FAIL);
      }
    }
 
    free(string_product_type);

/************/

    if (!strcmp(rp->type,"CPX"))
        sprintf(job_name,"REV %s %.05dX%.05d", platform, asp_rev_id,
		asp_job_id);
    else
        sprintf(job_name,"REV %s %.05d%.1s%.05d", platform, asp_rev_id,
		&rp->type[0], asp_job_id);
    free(platform);
    asp_msg(LOG_DEBUG,"%s", job_name); 

    crp->product_id = ODLGetString(rqst, PRODUCT_ID, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", PRODUCT_ID); 
	return(FAIL);
    }

    if (!strcmp(rp->type,"NUC")) {
	ucs_product = 1;	/* set flag of uncompensated */
	strcpy(rp->type,"RPR");
    }

    if (!strcmp(rp->type,"QLK")) { 
	quicklook = 1;	/* set flag of quicklook*/
	strcpy(rp->type,"RPR");
    }

    crp->version= ODLGetInt(rqst, FILE_VERSION, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", FILE_VERSION);
	return(FAIL);
    }
    if (vbose) printf("file_version = %d\n", crp->version);

/* ----------- Combine Frame ID, Subfram ID, and Res-Type ----------- */

    memset(rp->site,SP,sizeof(rp->site));

    crp->pix_spacing = (float) ODLGetDouble( rqst, PIXEL_SPACING, &status );
    for ( i = 0; (cp_pix_spacing[i] != NULL) &&
			(crp->pix_spacing != cp_pix_spacing[i]); ++i );
    /*if ( cp_pix_spacing[i] == NULL ){ 
	asp_msg(LOG_DEBUG,"invalid %s", PIXEL_SPACING);
	return(FAIL);
    }*/
    if ( strcmp(rp->type,asp_product_type[RPR]) && (!GetResType()) ){
	asp_msg(LOG_DEBUG,"ignore %s", PIXEL_SPACING);
	crp->pix_spacing = cp_pix_spacing[FULL_RES];
    }
    if (!strcmp(rp->type,asp_product_type[RPR])){
      if ( cp_pix_spacing[i] == NULL ){ 
	asp_msg(LOG_DEBUG,"invalid %s", PIXEL_SPACING);
	return(FAIL);
      }
    }
    crp->frame_ID = ODLGetInt(rqst, STD_FRAME_ID, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", STD_FRAME_ID);
	return(FAIL);
    }
    crp->subframe_ID = ODLGetInt(rqst, CPX_FRAME_ID, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", CPX_FRAME_ID);
	return(FAIL);
    } 
    rp->lat = rp->lon = rp->targ_rg = 0.;
    rp->ave_hght = (float) ODLGetDouble( rqst, AVG_TERRAIN_HT,
					&status );
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", AVG_TERRAIN_HT);
	return(FAIL);
    }
    rp->proc_gain = ODLGetInt(rqst, PROCESSING_GAIN, &status);
    if ( status == FAIL ) {
	asp_msg(LOG_DEBUG,"can't get %s", PROCESSING_GAIN);
	return(FAIL);
    }
    if (vbose) printf("Got rp->ave_hght %g rp->proc_gain %d\n",
		rp->ave_hght, rp->proc_gain);
    string = ODLGetString(rqst, DESKEW_FLAG, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", DESKEW_FLAG);
	free(string);
	return(FAIL);
    } else {
	strcpy( rp->deskew, string );
	if ( !strcmp(rp->deskew,"NO") ) strcpy(rp->deskew,"NOT");
    }
    free(string);

    string = ODLGetString(rqst, PROJECTION, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", PROJECTION);
        free(string);
	return(FAIL);
    } else {
	for ( i = 0; (cp_projection[i] != NULL) &&
			strcmp(string,cp_projection[i]); ++i );
	if ( cp_projection[i] != NULL ) 
			strcpy( rp->gnd_slnt_rg, asp_projection[i] );
	else {
		asp_msg(LOG_DEBUG,"invalid %s", PROJECTION);
        	free(string);
		return(FAIL);
	}
        free(string);
    }
    if (vbose) printf("Got rp->deskew %s rp->gnd_slnt_rg %s\n",
			rp->deskew, rp->gnd_slnt_rg);

/* ----------- Get input/output filenames */

/* ASP-CP 2.1, 11/21/96 */

    crp->cal_file = ODLGetString(rqst, CALPARMS_FILE_PRIMARY, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", CALPARMS_FILE_PRIMARY);
	return(FAIL);
    }
/*
    crp->cal_file = ODLGetString(rqst, PVS_CALPARMS_FILE, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", PVS_CALPARMS_FILE);
	return(FAIL);
    }
*/
    if (vbose)
	printf("cal_file %s\n", crp->cal_file);
    crp->image_file = ODLGetString(rqst, IMAGE_FILE, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", IMAGE_FILE);
	return(FAIL);
    }
    if (vbose)
	printf("image_file %s\n", crp->image_file);
    crp->ceos_ldr = ODLGetString(rqst, CEOS_LEADER_FILE, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", CEOS_LEADER_FILE);
	return(FAIL);
    }
/* Comment out: extract version # from the sps_frame_rqst CV 3/6/96 
    } else {
	strncpy(stmp,&crp->ceos_ldr[13],3);	
	crp->version = atoi(stmp);
	if ( (crp->version <= 0) || (crp->version >= 1000) ){
		if (vbose) printf("version not within range of 3 digits\n");
		crp->version = 1;
	}
    }
*/
    if (vbose)
	printf("ceos_ldr %s\n", crp->ceos_ldr);
    crp->pmf_file = ODLGetString(rqst, PMF_FILE, &status);
    if ( status == FAIL ){
	asp_msg(LOG_DEBUG,"can't get %s", PMF_FILE);
	return(FAIL);
    }
    if (vbose)
	printf("pmf_file %s\n", crp->pmf_file);

    if (!strcmp(rp->type,asp_product_type[RPR])) { 
	if ( GetResType() ) {
	   if (!strncmp(&crp->product_id[12],"U",1)){
		printf("TEST\n");
                sprintf(rp->site, "%3.3d%0.1s1U%3.3d", crp->frame_ID,
			rp->gnd_slnt_rg, crp->version ); 
	   } else if (ucs_product) { 
                sprintf(rp->site, "%3.3d%0.1s1N%3.3d", crp->frame_ID,
			rp->gnd_slnt_rg, crp->version ); 
           } else if (quicklook) {
		sprintf(rp->site, "%3.3d%0.1s1Q%3.3d", crp->frame_ID,
			rp->gnd_slnt_rg, crp->version ); 
           } else {
		sprintf(rp->site, "%3.3d%0.1s1S%3.3d", crp->frame_ID,
			rp->gnd_slnt_rg, crp->version ); 
	   }
	} else {
	   if (!strncmp(&crp->product_id[12],"U",1)){
		printf("TEST\n");
                sprintf(rp->site, "%3.3d%0.1s4U%3.3d", crp->frame_ID,
			rp->gnd_slnt_rg, crp->version ); 
	   } else if (ucs_product) { 
                sprintf(rp->site, "%3.3d%0.1s4N%3.3d", crp->frame_ID,
			rp->gnd_slnt_rg, crp->version ); 
           } else if (quicklook) {
		sprintf(rp->site, "%3.3d%0.1s4Q%3.3d", crp->frame_ID,
			rp->gnd_slnt_rg, crp->version ); 
           } else {
		sprintf(rp->site, "%3.3d%0.1s4S%3.3d", crp->frame_ID,
			rp->gnd_slnt_rg, crp->version ); 
	   }
	}
    } else if ( !strcmp(rp->type,asp_product_type[CSD]) ){
        sprintf(rp->site, "%3.3d%0.1s0C%3.3d", crp->frame_ID,  
			rp->gnd_slnt_rg, crp->version ); 
    } else if ( !strcmp(rp->type,asp_product_type[CPX]) ){
        sprintf(rp->site, "%3.3d%1.1d0X%3.3d",crp->frame_ID,crp->subframe_ID,
			crp->version ); 
    }

    return(GetODLRqst(rqst,crp,rp,sv,tcf));

} /* end of GetFrameRqst */

int GetODLGMT( rqst, name, gmt )
	ODL rqst;
	char *name;
	GMT_PTR gmt;
{
	char *string;
	int status, is, fs;

	if ( (string = ODLToStr(rqst,name)) == NULL ) {
		asp_msg(LOG_DEBUG,"can't get %s", name);
		free(string);
		return(FAIL);
    	} else {
		if (vbose) printf("string = %s, ", string);
		strcpy(string,strtok(string,"="));
		strcpy(string,strtok(NULL,"\0"));
		if (vbose) printf("string = %s\n", string);
		get_gmt( string, gmt );
		free(string);
    	}
    	if (vbose) {
		printf("%s:",name);
 	        is = (int)(gmt->second);
        	fs = (int)((gmt->second - (float)(is))*1000.0+0.5);
		printf("%04.4d:%03.3d:%02.2d:%02.2d:%02.2d.%03.3d\n",gmt->yr,
			gmt->day, gmt->hr, gmt->min, is,fs);
	}
	return(PASS);
} /* end of GetODLGMT */

ODL GetODLFile( filename )
	char *filename;
{
	char err[256];
	ODL odl;

	asp_msg(LOG_DEBUG,"GetODLFile: reading %s",filename); 

	if ((odl = ODLparse(filename, 0, err)) == NULL)
	    printf("can't parse %s, %s. ODL is null\n", filename, err);

	return( odl );

} /* end of GetODLFile */

int GetODLScanFile ( frame_file, frame_ID, sp_out, pre_out, rp_out )
	char *frame_file;
	int frame_ID;
	TAPE_SEG_PTR *sp_out;
	PREAM_FILE_PTR *pre_out;
	RQST_PTR *rp_out;	
{
	ODL	odl;
	int	status, i, j, m, n;
	int	seg_cnt, pp_cnt, frame_count, cur_frame_ID;
	char	**msg;
	TAPE_SEG_PTR sp, GetODLSegment();
	PREAM_FILE_PTR pre, GetODLPreamble();
	ODL *odl1,*segmt, *frame, *agcchg, *agcfmt;
	char string[132];
	int st_fmt_next;
        int block, offset, norm_inv, itime;
        int ami_fmt, pri, format, fst_len;
	double prf;
	int start_fmt, end_fmt;

	if ( (odl = GetODLFile( frame_file )) == NULL ||
	     (segmt = (ODL*) Val(Lookup(odl, SCAN_BODY))) == NULL ) {
	   if ( (segmt = (ODL*) Val(Lookup(odl, "BODY"))) == NULL ) {
	      asp_msg(LOG_DEBUG,"GetODLScanFile: can't read %s", frame_file);
	      ODLFree(odl); return(FAIL);
	   }
	}
	src= ODLGetString(odl,SCAN_SOURCE,&status);
	if (status == FAIL) {
           asp_msg(LOG_DEBUG,"GetODLScanFile: can't read %s",SCAN_SOURCE);
           ODLFree(odl); return(FAIL);
	}
	for (seg_cnt = 0,i = 0; segmt[i] != NULL; ++i)
	{  
	   if (!strcasecmp(Name(segmt[i]), PRE_CAL1_POW)) {
		cal.pre_cal1_pow = 
			(float) ODLGetDouble( segmt[i], PRE_CAL1_POW, &status );
		continue;
	   } else if (!strcasecmp(Name(segmt[i]), PRE_CAL2_POW)) {
		cal.pre_cal2_pow = 
			(float) ODLGetDouble( segmt[i], PRE_CAL2_POW, &status );
		continue;
	   } else if (!strcasecmp(Name(segmt[i]), POST_CAL1_POW)) {
		cal.post_cal1_pow = 
			(float)ODLGetDouble(segmt[i], POST_CAL1_POW, &status);
		continue;
	   } else if (!strcasecmp(Name(segmt[i]), POST_CAL2_POW)) {
		cal.post_cal2_pow = 
			(float)ODLGetDouble(segmt[i], POST_CAL2_POW, &status);
		continue;
	   } else if (strcasecmp(Name(segmt[i]), SEGMENT))
		continue;

	    /* loop for SEGMENT */

	   ++seg_cnt;
	   frame_count = ODLGetInt(segmt[i], FRAME_COUNT, &status); 
	   if ( status == FAIL ){
		asp_msg(LOG_DEBUG,
			"GetODLScanFile: can't not read %s in segment %d",
			FRAME_COUNT, seg_cnt); 
		ODLFree(odl); return(FAIL);
	   }
	   if (vbose) 
		printf("Segment %d has %d frames\n", seg_cnt, frame_count);
	   
	   if ( (pre = GetODLPreamble(segmt[i])) == NULL ){
		asp_msg(LOG_DEBUG,"GetODLScanFile: GetODLPreamble failed");
		ODLFree(odl); return(FAIL);
	   }
	   if ( (sp = GetODLSegment(segmt[i],pre)) == NULL ){
		asp_msg(LOG_DEBUG,"GetODLScanFile: GetODLSegment failed");
		ODLFree(odl); return(FAIL);
	   }
	   sp->pre = pre;

	   if ((frame = (ODL*) Val(segmt[i])) == NULL){
		asp_msg(LOG_DEBUG,
			"GetODLScanFile: can't find any %s\n", IMAGE_FRAME);
		ODLFree(odl); ODLFree(msg); return(FAIL);
	   }
	   for (pp_cnt = 0, j = 0; frame[j] != NULL; ++j)
	   {
	        if (strcasecmp(Name(frame[j]), IMAGE_FRAME))
		    continue;

	        /* loop for IMAGE_FRAME */
		odl1 = frame[0];
		cur_frame_ID = ODLGetInt(frame[j], FRAME_ID, &status);
		if ( vbose ) printf("GetODLScanFile: @frame id %d\n",
					cur_frame_ID);
		if ( status == FAIL ){
		   asp_msg(LOG_DEBUG,"GetODLScanFile: can't read %s", FRAME_ID);
		   ODLFree(odl); return(FAIL);
		}
		if ( !pp_cnt ) {
		   if ( FixSegInfo(frame[j],sp,pre) == FAIL ){
			asp_msg(LOG_DEBUG,"GetODLScanFile: FixSegInfo failed");	
			ODLFree(odl); return(FAIL);
		   }
		}
		if ( cur_frame_ID == frame_ID ) {
		  if (vbose) 
			printf("GetODLScanFile: found frame_ID %d\n", frame_ID);
        	  if (chkrds) {
		     if ((!strcmp(src,"RDS")) || (strcmp(media_type,"DCRSI"))) {
                        /* read the start_format & end_format from RDS scan
                            result file CV 7/24/97  */
                        start_fmt = ODLGetInt(frame[j], FRAME_ST_FMT, &status);
                        if ( status == FAIL ){
                             asp_msg(LOG_DEBUG,"GetODLScanFile: can't read %s",
                                FRAME_ST_FMT);
                             ODLFree(odl); return(FAIL);
                        }
           	        if (!strcmp(src,"RDS")) {
                          end_fmt = ODLGetInt(frame[j], FRAME_END_FMT, &status);
                          if ( status == FAIL ){
                             asp_msg(LOG_DEBUG,"GetODLScanFile:can't read %s",
                                FRAME_END_FMT);
                             ODLFree(odl); return(FAIL);
			  }
                          rds_center_format = (int)((start_fmt+end_fmt)/2);
                          if (vbose) {
			     printf("start_fmt=%d,end_fmt=%d,center=% d\n",
                                start_fmt, end_fmt, rds_center_format);
                          }
                        } /* RDS */
		        if ( VerifySync(frame[j],sp,pre) == FAIL ){
			   asp_msg(LOG_DEBUG,"GetODLScanFile: VerifySync fail\n");	
			   ODLFree(odl); return(FAIL);
		        }
		     }
        	  }
		  st_fmt = ODLGetInt(frame[j], FRAME_ST_FMT, &status);
		  if (vbose) printf("st_fmt=%d \n", st_fmt);
		  if ( status == FAIL ){
		      asp_msg(LOG_DEBUG,"GetODLScanFile: can't read %s", 
				FRAME_ST_FMT);
		      ODLFree(odl); return(FAIL);
		  }
		  max_agc = 0;
		  if (!strcmp(src,"RDS") && (Cur_Rqst->take_id[0] == 'R')) {
		  /* extract the right agc value aroung the frame processing
		  and write to pp_region file CV 4/9/97 */
		      if (frame[j+1] != NULL) {
		        st_fmt_next= ODLGetInt(frame[j+1],FRAME_ST_FMT,&status);
		        if ( status == FAIL )
		            st_fmt_next = st_fmt + 8000;
			printf("st_fmt_next=%d, ", st_fmt_next);
		      }		
	              if ((agcchg = (ODL*) Val(Lookup(segmt[i], AGC_VALUE_CHANGE))) == NULL ){
	                  asp_msg(LOG_DEBUG, "GetODLScanFile: can't read %s",
				AGC_VALUE_CHANGE);
                	  return(NULL);
        	      }
        	      if ((agcfmt = (ODL*) Val(Lookup(segmt[i], AGC_FORMAT_CHANGE))) == NULL){
                	  asp_msg(LOG_DEBUG, "GetODLScanFile: can't read %s",
			 	AGC_FORMAT_CHANGE);
                	  return(NULL);
        	      }
        	      /*for ( m=0,n=0; (agcchg[m] != NULL) && (n<80); ++m ){*/
        	      for ( m=0,n=0; (agcchg[m] != NULL) && ((*(int*) Val(agcfmt[m])) < st_fmt_next); ++m ){
			  if ((*(int*) Val(agcfmt[m])) > st_fmt) { 
			    if (n < 80) {
                	       sp->agc[n] = *(int*) Val(agcchg[m-1]);
                	       sp->afmt[n] = *(int*) Val(agcfmt[m-1]);
			       n++;
			    }
			    if ((*(int*) Val(agcchg[m-1])) > max_agc)
                	       max_agc = *(int*) Val(agcchg[m-1]);
			  }
        	      }
        	      sp->agc_count = n;
                      if (sp->agc_count == 0) {
                          printf("No agcfmt is greater than st_fmt %d\n",m);
                          sp->agc[n] = *(int*) Val(agcchg[m-1]);
                          sp->afmt[n] = *(int*) Val(agcfmt[m-1]);
                          sp->agc_count += 1;
                          max_agc = sp->agc[n];
                      }
		   }
		   else if (!strcmp(src,"ASP") && (Cur_Rqst->take_id[0] == 'R')) {
		 	for (m=0; m<sp->agc_count; m++) {
			   if (sp->agc[m] > max_agc)
			   	max_agc = sp->agc[m];
			}
		   }
        	   if (vbose) printf("m=%d,agc_count=%d, PEAK AGC=%d\n",
			m,sp->agc_count,max_agc);
		   if (max_agc > 32)
			printf("WARNING: AGC IS GREATER THAN 32 =%d\n",max_agc);

		   if ( GetODLFrame(frame[j],rp_out) == FAIL ){
			asp_msg(LOG_DEBUG,"GetODLScanFile: GetODLFrame failed");
			ODLFree(odl); return(FAIL);
		   }
		   *sp_out = sp; *pre_out = pre;
		   save_seg = (char *) ODLToStr(segmt[i],NULL);
		   save_frame = (char *) ODLToStr(frame[j],NULL);
		   ODLFree(odl); return(PASS);
		} 
		pp_cnt++;
	   } /* end of frame loop */
	   if ( pp_cnt != frame_count ){
	      asp_msg(LOG_DEBUG,"WARNING: segment %d has inconsistent number of frames");
	      asp_msg(LOG_DEBUG,"%s %d counter %d",FRAME_COUNT,frame_count,
				pp_cnt);
	   }
	   free(sp); free(pre);
	} /* end of segment loop */
	asp_msg(LOG_ERR,asp_messages[UNFOUND_FRAME_ID]);
	return(FAIL);
} /* end of GetODLScanFile */

TAPE_SEG_PTR GetODLSegment(odl)
	ODL odl;
{
	TAPE_SEG_PTR sp;
	int status, i, count;
	ODL	*window, *format, *agcchg, *agcfmt;
	char *string;

	if ( (sp = (TAPE_SEG_PTR) malloc(sizeof(TAPE_SEG))) == NULL ){
		asp_msg(LOG_DEBUG,"GetODLSegment: can't allocate TAPE_SEG\n");
		return(NULL);
	}
	sp->ppb_list = NULL;
	sp->nxt_seg = NULL;
	sp->post = NULL;
	memset(sp,0,sizeof(TAPE_SEG));

	if ( Cur_Rqst->take_id[0] == 'R' ){
	        string = ODLGetString(odl, SEG_MODE, &status);
		for ( i = 1; (R_mode[i] != NULL) && 
				strcmp(string,R_mode[i]); ++i );
		free(string);
		if ( strcmp(R_mode[i], cp_rp->mode)) {
		   asp_msg(LOG_DEBUG,"MODE in job_file & scan_result are MISMATCH\n"); 
		   return(NULL);
		}
		if ( R_mode[i] != NULL ) sp->aux.beam_seq = i;
		else {
			asp_msg(LOG_DEBUG,"unsupported %s", SEG_MODE);
			return(NULL);
		}
	} else	sp->aux.beam_seq = 0;
	if ( vbose ) printf("RSAT MODE %d\n",sp->aux.beam_seq);

	if (GetODLGMT(odl, SEG_START_TIME, &sp->start) == FAIL ) 
			return(NULL);
	if (GetODLGMT(odl, SEG_END_TIME, &sp->end) == FAIL ) 
			return(NULL); 

	sp->blk_start = ODLGetInt(odl, SEG_START_TAPE_BLOCK, &status);
	if ( status == FAIL ){
		asp_msg(LOG_DEBUG,"GetODLSegment: can't read %s", 
				SEG_START_TAPE_BLOCK);
		return(NULL);
	}
	sp->blk_end = ODLGetInt(odl, SEG_END_TAPE_BLOCK, &status);
	if ( status == FAIL ){
		asp_msg(LOG_DEBUG,"GetODLSegment: can't read %s", 
				SEG_END_TAPE_BLOCK);
		return(NULL);
	}
	if ((format = (ODL*) Val(Lookup(odl, FORMAT_CHANGE))) == NULL ){
		asp_msg(LOG_DEBUG,
			"GetODLSegment: can't read %s", FORMAT_CHANGE);
		return(NULL);
	}
	if ((window = (ODL*) Val(Lookup(odl, WINDOW_POSN_CHANGE))) == NULL ){
		asp_msg(LOG_DEBUG,
			"GetODLSegment: can't read %s", WINDOW_POSN_CHANGE);
		return(NULL);
	}
	for ( i = 0; (window[i] != NULL) && (i < 20); ++i ){ 
		sp->wfmt[i] = *(int*) Val(format[i]);
		sp->wdwp[i] = *(int*) Val(window[i]);
	}
	sp->win_count = i;
	if (vbose) printf("win_count %d\n",sp->win_count);

	if ((agcchg = (ODL*) Val(Lookup(odl, AGC_VALUE_CHANGE))) == NULL ){
		asp_msg(LOG_DEBUG,
			"GetODLSegment: can't read %s", AGC_VALUE_CHANGE);
		return(NULL);
	}
	if ((agcfmt = (ODL*) Val(Lookup(odl, AGC_FORMAT_CHANGE))) == NULL ){
		asp_msg(LOG_DEBUG,
			"GetODLSegment: can't read %s", AGC_FORMAT_CHANGE);
		return(NULL);
	}
	for ( i = 0; (agcchg[i] != NULL) && (i < 80); ++i ){ 
		sp->agc[i] = *(int*) Val(agcchg[i]);
		sp->afmt[i] = *(int*) Val(agcfmt[i]);
	}
	sp->agc_count = i;
	if (vbose) printf("agc_count %d\n",sp->agc_count);
	sp->polarity = cp_rp->activity;

	return(sp);

} /* end of GetODLSegment */

PREAM_FILE_PTR GetODLPreamble(odl)
	ODL odl;
{
	PREAM_FILE_PTR sp;
	int status;

	if ( (sp = (PREAM_FILE_PTR) malloc(sizeof(PREAM_FILE))) == NULL ){
		asp_msg(LOG_DEBUG,"GetODLPreamble: can't allocate PREAM_FILE");
		return(NULL);
	}
	memset(sp,0,sizeof(PREAM_FILE));
	sp->nxt_pre = NULL;
	strncpy( sp->sat, Cur_Rqst->take_id, 2 ); 
	strcpy( sp->takeid, Cur_Rqst->take_id );
	strcpy( sp->dcrs_id, Cur_Rqst->tape_id );
	strcpy( sp->sen_mode, "OGRC" );
	sp->prf = ODLGetDouble( odl, SCAN_PRF, &status );
	sp->polarity = cp_rp->activity;
	return(sp);

} /* end of GetODLPreamble */

int VerifySync (odl, sp, pre)
	ODL odl;
	TAPE_SEG_PTR sp;
	PREAM_FILE_PTR pre;
{

	int found_sync, rec1, nblks;
        int block, offset, norm_inv, itime;
        int ami_fmt, pri, format, fst_len;
        int start_fmt, end_fmt,no_block;
        int st_blk, st_bit_offset, status, end_blk,st_fmt;
        double est_bits, est_bits2;
	double total_bit, total_bytes, m_bytes;
        double prf;
	int TAPE_BITS, TAPE_BLOCK;
	int cur_frame_id;
	int count = 0;

	found_sync = 0;
        if (!strcmp( media_type, "DCRSI" ) ) {
             TAPE_BITS = 4356*8;
             TAPE_BLOCK = 4356;
             nblks = 60;
        } else {
             TAPE_BITS = 144304*8;
             TAPE_BLOCK = 144304;
             nblks = 2;
        }

        if (!strcmp( data_direction, "REVERSE" ) )
            block = sp->blk_end - 20; /* CV 9/23/97 */
        else
            block = sp->blk_start;
                                                       

	st_blk = ODLGetInt(odl, FRAME_ST_BLK, &status);
	if ( status == FAIL ){
		asp_msg(LOG_DEBUG,"can't read %s", FRAME_ST_BLK);
		return(FAIL);
	}
	st_bit_offset = ODLGetInt(odl, FRAME_ST_BIT, &status);
	if ( status == FAIL ){
		asp_msg(LOG_DEBUG,"can't read %s", FRAME_ST_BIT);
		return(FAIL);
	}
	st_fmt = ODLGetInt(odl, FRAME_ST_FMT, &status);
	if ( status == FAIL ){
		asp_msg(LOG_DEBUG,"can't read %s", FRAME_ST_FMT);
		return(FAIL);
	}
	end_blk = ODLGetInt(odl, FRAME_END_BLK, &status);
	if ( status == FAIL ){
		asp_msg(LOG_DEBUG,"can't read %s", FRAME_ST_BLK);
		return(FAIL);
	}
        cur_frame_id = ODLGetInt(odl, FRAME_ID, &status);
        if ( vbose ) printf("VerifySync: frame id %d\n", cur_frame_id);
        if ( status == FAIL ){
            asp_msg(LOG_DEBUG,"VerifySync: can't read %s", FRAME_ID);
            ODLFree(odl); return(FAIL);
        }
printf("VUU st_blk=%d,st_bit_offset=%d,st_fmt=%d,end_blk=%d\n",
	  st_blk, st_bit_offset, st_fmt, end_blk);
        if ( Cur_Rqst->take_id[0] == 'E' ){
          if ((!strcmp(data_direction,"REVERSE")) &&
                (!strcmp(chk_media_id,"YES"))) {
              if (e1_srch_sync(3,2,nblks,&sp->blk_start,&offset,
                  &norm_inv,&format,&ami_fmt,&pri) != PASS) {
                 if (Write_CP_Error_Status == -2)
                    return(FAIL);
              } else
                 strcpy(chk_media_id,"NO");
           }
           if (e1_srch_sync(3,2,nblks,&block,&offset,
                  &norm_inv,&format,&ami_fmt,&pri) != PASS) {
              if (Write_CP_Error_Status == -2) {
                 return(FAIL); /* tape_media_id is not matched */
              } else {
                 strcpy(chk_media_id,"NO");
                 while (count < 2) {
                    if( strcmp( media_type, "DCRSI" ) ){ /* SONY */
                       sony_stop(IN);
                       block += 15;
                       if (block > sp->blk_end) block -= 15;
                    } else { 
                       dc_stop(IN);
                       block += 500;
		    }
                    if (e1_srch_sync(3,2,nblks,&block,&offset,
                       &norm_inv,&format,&ami_fmt,&pri) != PASS){
                       count += 1;
                       if (count == 2 ) {
                          printf("RDS SYNC_OFFSET IS OFF TOO MUCH\n ");
                          return (FAIL); /* CV 9/23/97, no srch frame sync */
                       }
                    }
                    else {
                       found_sync = 1;
                       break;
                    }
                 }
              }
              if( strcmp( media_type, "DCRSI" ) ) /* SONY */
                 sony_stop(IN);
              else
                 dc_stop(IN);
           } else found_sync = 1;
        } else if ( Cur_Rqst->take_id[0] == 'J' ){
           if ((!strcmp(data_direction,"REVERSE")) &&
                (!strcmp(chk_media_id,"YES"))) {
              if (j1_srch_sync(3,2,nblks,&sp->blk_start,&offset,
                  &norm_inv,&format,&prf,&fst_len) != PASS) {
                 if (Write_CP_Error_Status == -2) 
                    return(FAIL); /* tape_media_id is not matched */
              } else 
                    strcpy(chk_media_id,"NO");
	   }		
           if (j1_srch_sync(3,2,nblks,&sp->blk_start,&offset,
                  &norm_inv,&format,&prf,&fst_len) != PASS) {
              if (Write_CP_Error_Status == -2) { 
                  return(FAIL); /* tape_media_id is not matched */
              } else {
                  strcpy(chk_media_id,"NO");
                  while (count < 2) {
                     if( strcmp( media_type, "DCRSI" ) ){ /* SONY */
                        sony_stop(IN);
                        block += 15;
                        if (block > sp->blk_end) block -= 15;
                     } else { 
                        dc_stop(IN);
                        block += 500;
		     }
                     if (j1_srch_sync(3,2,nblks,&block,&offset,
                          &norm_inv,&format,&prf,&fst_len) != PASS){
                        count += 1;
                        if (count == 2 ) {
                            printf("RDS SYNC_OFFSET IS OFF TOO MUCH\n ");
                            return (FAIL); /* CV 9/23/97, no srch frame sync */
                        }
                      } else {
                        found_sync = 1;
                        break;
                      }
                   }
              }
              if( strcmp( media_type, "DCRSI" ) ) /* SONY */
                 sony_stop(IN);
              else
                 dc_stop(IN);
           } else found_sync = 1;
        } else if ( Cur_Rqst->take_id[0] == 'R' ){
          if ((!strcmp(data_direction,"REVERSE")) &&
                (!strcmp(chk_media_id,"YES"))) {
              if (r1_srch_sync(3,2,nblks,&sp->blk_start,&offset,
                  &norm_inv,&itime,&prf) != PASS) {
                 if (Write_CP_Error_Status == -2)
                    return(FAIL);
              } else
                 strcpy(chk_media_id,"NO");
           }
           if (r1_srch_sync(3,2,nblks,&block,&offset,
                  &norm_inv,&itime,&prf) != PASS) {
              if (Write_CP_Error_Status == -2) {
                 return(FAIL); /* tape_media_id is not matched */
              } else {
                 strcpy(chk_media_id,"NO");
                 while (count < 2) {
                    if( strcmp( media_type, "DCRSI" ) ){ /* SONY */
                       sony_stop(IN);
                       block += 15;
                       if (block > sp->blk_end) block -= 15;
                    } else { 
                       dc_stop(IN);
                       block += 500;
		    }
                    if (r1_srch_sync(3,2,nblks,&block,&offset,
                       &norm_inv,&itime,&prf) != PASS){
                       count += 1;
                       if (count == 2 ) {
                          printf("RDS SYNC_OFFSET IS OFF TOO MUCH\n ");
                          return (FAIL); /* CV 9/23/97, no srch frame sync */
                       }
                    }
                    else {
                       found_sync = 1;
                       break;
                    }
                 }
              }
              if( strcmp( media_type, "DCRSI" ) ) /* SONY */
                 sony_stop(IN);
              else
                 dc_stop(IN);
           } else found_sync = 1;
           if (sp->polarity != norm_inv) {
              if (vbose) {
		printf("POLARITY =%d\n",norm_inv);
		printf("ACTIVITY_ID in job file is INCORRECT\n");	
	      }
	      return(FAIL);
	   }
        }
	if (found_sync) {
           printf("found sync at block = %d offset %d\n", block, offset);
           if ( Cur_Rqst->take_id[0] == 'R' ){
             if ((((abs(offset-sp->bit_off))/8)%323) != 0)
                printf("RDS offset is off, ");
	   }
           if ( use_rds_scan ){
              printf("SEGMENT : calculate: %d %d, srch: %d %d\n",
                       sp->blk_start, sp->bit_off, block, offset);
              if (!strcmp( data_direction, "REVERSE" ) )
                   sp->blk_end = block;
              else
                   sp->blk_start = block;
              if ((!strcmp(host,"newasp")) && (strcmp(media_type,"DCRSI")))
                 sp->bit_off = offset + 8;
              else
                 sp->bit_off = offset;
	   }
	   return(PASS);
        }
}

int FixSegInfo( odl, sp, pre )
	ODL odl;
	TAPE_SEG_PTR sp;
	PREAM_FILE_PTR pre;
{
	double deltat, get_gmt_diff();
	int st_blk, st_bit_offset, status, end_blk;
	double est_bits, est_bits2;
	GMT st_gmt;
	int TAPE_BITS;
printf("VUU pass here\n");
	if ( Cur_Rqst->take_id[0] == 'E' ){
		sp->fmt_id = sp->wfmt[0];
		sp->fmt_start = -1;
	} else {
		sp->fmt_start = sp->wfmt[0];
	}
	deltat = get_gmt_diff( &sp->end, &sp->start );
printf("VUU pass here 1\n");
	if ( Cur_Rqst->take_id[0] == 'E' ){
	   sp->end_id = sp->fmt_id + deltat * pre->prf;
           pre->time_ref_fmt = sp->fmt_id;
	} else {
	   sp->fmt_end = sp->fmt_start + deltat * pre->prf;
           pre->time_ref_fmt = sp->fmt_start;
	}
        pre->time_ref_gmt = sp->start;

	st_blk = ODLGetInt(odl, FRAME_ST_BLK, &status);
	if ( status == FAIL ){
		asp_msg(LOG_DEBUG,"can't read %s", FRAME_ST_BLK);
		return(FAIL);
	}
        end_blk = ODLGetInt(odl, FRAME_END_BLK, &status);
        if ( status == FAIL ){
                asp_msg(LOG_DEBUG,"can't read %s", FRAME_ST_BLK);
                return(FAIL);
        }
	st_bit_offset = ODLGetInt(odl, FRAME_ST_BIT, &status);
	if ( status == FAIL ){
		asp_msg(LOG_DEBUG,"can't read %s", FRAME_ST_BIT);
		return(FAIL);
	}
	st_fmt = ODLGetInt(odl, FRAME_ST_FMT, &status);
	if ( status == FAIL ){
		asp_msg(LOG_DEBUG,"can't read %s", FRAME_ST_FMT);
		return(FAIL);
	}
        if (GetODLGMT(odl, "START_TIME", &st_gmt) == FAIL ){ 
		asp_msg(LOG_DEBUG,"can't read %s", "START_TIME");
		return(FAIL);
	}
/*SONY*/
        if (!strcmp( media_type, "DCRSI" ) ) 
             TAPE_BITS = 4356*8;
        else 
             TAPE_BITS = 144304*8;
printf("VUU pass here 2\n");

first_st_blk = st_blk;
first_st_bit = st_bit_offset;	
#define E_FMT_BITS (29*256*8)
#define R_FRAME_BITS (323*8)

        if (!strcmp( data_direction, "REVERSE" ) )
            est_bits = (double) (sp->blk_end - end_blk) * TAPE_BITS +
                                st_bit_offset;
        else
	    est_bits = (double) (st_blk - sp->blk_start) * TAPE_BITS +
				st_bit_offset; 
	if ( Cur_Rqst->take_id[0] == 'E' ){
	   est_bits2 = (st_fmt - sp->fmt_id) * data_rate[sp->polarity] / 
			pre->prf;
	   sp->bit_off = est_bits - nint(est_bits2/E_FMT_BITS)*E_FMT_BITS; 
	   while ( sp->bit_off <= 0 ) sp->bit_off += E_FMT_BITS;
           while ( sp->bit_off > E_FMT_BITS ) sp->bit_off -=
E_FMT_BITS;
	} else if ( Cur_Rqst->take_id[0] == 'R' ) {
	   est_bits2 = (st_fmt - sp->fmt_start) * (data_rate[sp->polarity] +
840)/ pre->prf;
	   sp->bit_off = est_bits - nint(est_bits2/R_FRAME_BITS)*R_FRAME_BITS; 
	   while ( sp->bit_off <= 0 ) sp->bit_off += R_FRAME_BITS;/*CV 4/14/97*/

/* fake consistent blk & fmt in pp_regions for J1_get_fmt_loc 02/13/97 */
	} else if ( Cur_Rqst->take_id[0] == 'J' ) { 
	   sp->fmt_start = st_fmt;
	   sp->bit_off = st_bit_offset;
	   sp->blk_start = st_blk;
           pre->time_ref_fmt = sp->fmt_start;
           pre->time_ref_gmt = st_gmt;
	   sp->start = st_gmt;
	}
printf("VUU pass here 3\n");

} /* end of FixSegInfo */
int GetODLFrame (odl, rp)
	ODL odl;
	RQST_PTR rp;	
{
	int status;
	double a_lat,b_lat,c_lat,d_lat,e_lat;
	double a_lon,b_lon,c_lon,d_lon,e_lon;
	float index;

	if (!strcmp(Cur_Rqst->type,"CPX")) {
	    c_lat = (float) ODLGetDouble(odl, CENTER_LAT, &status);
	    if ( status == FAIL ){
		asp_msg(LOG_DEBUG,"GetODLFrame: can't read %s", CENTER_LAT);
		return(FAIL);
	    }
	    c_lon = (float) ODLGetDouble(odl, CENTER_LON, &status);
	    if ( status == FAIL ){
		asp_msg(LOG_DEBUG,"GetODLFrame: can't read %s", CENTER_LON);
		return(FAIL);
	    }
	    d_lat = (float) ODLGetDouble(odl, NEAR_START_LAT, &status);
	    if ( status == FAIL ){
	      asp_msg(LOG_DEBUG,"GetODLFrame: can't read %s", NEAR_START_LAT);
	      return(FAIL);
	    }
	    d_lon = (float) ODLGetDouble(odl, NEAR_START_LON, &status);
	    if ( status == FAIL ){
	      asp_msg(LOG_DEBUG,"GetODLFrame: can't read %s", NEAR_START_LON);
	      return(FAIL);
	    }
	    if (d_lon < 0 ) d_lon += 360;
	    a_lat = (float) ODLGetDouble(odl, NEAR_END_LAT, &status);
	    if ( status == FAIL ){
	      asp_msg(LOG_DEBUG,"GetODLFrame: can't read %s", NEAR_END_LAT);
	      return(FAIL);
	    }
	    a_lon = (float) ODLGetDouble(odl, NEAR_END_LON, &status);
	    if ( status == FAIL ){
	      asp_msg(LOG_DEBUG,"GetODLFrame: can't read %s", NEAR_END_LON);
	      return(FAIL);
	    }
	    if (a_lon < 0 ) a_lon += 360;
	    e_lat = (float) ODLGetDouble(odl, FAR_START_LAT, &status);
	    if ( status == FAIL ){
	      asp_msg(LOG_DEBUG,"GetODLFrame: can't read %s", FAR_START_LAT);
	      return(FAIL);
	    }
	    e_lon = (float) ODLGetDouble(odl, FAR_START_LON, &status);
	    if ( status == FAIL ){
	      asp_msg(LOG_DEBUG,"GetODLFrame: can't read %s", FAR_START_LON);
	      return(FAIL);
	    }
	    if (e_lon < 0 ) e_lon += 360;
	    b_lat = (float) ODLGetDouble(odl, FAR_END_LAT, &status);
	    if ( status == FAIL ){
	      asp_msg(LOG_DEBUG,"GetODLFrame: can't read %s", FAR_END_LAT);
	      return(FAIL);
	    }
	    b_lon = (float) ODLGetDouble(odl, FAR_END_LON, &status);
	    if ( status == FAIL ){
	      asp_msg(LOG_DEBUG,"GetODLFrame: can't read %s", FAR_END_LON);
	      return(FAIL);
	    }
	    if (b_lon < 0 ) b_lon += 360;
	    printf("subframe_ID=%d; ",cp_rp->subframe_ID);
/*    	    switch (cp_rp->subframe_ID) {
	        case 1: 
			rp->lat = a_lat + ((c_lat - a_lat)/2);
			rp->lon = a_lon + ((b_lon - a_lon)/8);
			break;
	        case 2: 
			rp->lat = a_lat + ((c_lat - a_lat)/2);
			rp->lon = a_lon + ((b_lon - a_lon)*3/8);
			break;
	        case 3: 
			rp->lat = a_lat + ((c_lat - a_lat)/2);
			rp->lon = a_lon + ((b_lon - a_lon)*5/8);
			break;
	        case 4: 
			rp->lat = a_lat + ((c_lat - a_lat)/2);
			rp->lon = a_lon + ((b_lon - a_lon)*7/8);
			break;
	        case 5: 
			rp->lat = c_lat;
			rp->lon = a_lon + ((b_lon - a_lon)/8);
			break;
	        case 6: 
			rp->lat = c_lat;
			rp->lon = a_lon + ((b_lon - a_lon)*3/8);
			break;
	        case 7: 
			rp->lat = c_lat;
			rp->lon = a_lon + ((b_lon - a_lon)*5/8);
			break;
	        case 8: 
			rp->lat = c_lat;
			rp->lon = a_lon + ((b_lon - a_lon)*7/8);
			break;
	    }
*/
	    if ( Cur_Rqst->take_id[0] == 'R' ) 
	     	index = 0.22;
	    if ( Cur_Rqst->take_id[0] == 'J' ) 
	    	index = 0.35;
	    if ( Cur_Rqst->take_id[0] == 'E' ) 
	    	index = 0.25;

	    if (cp_rp->subframe_ID < 5) {
		rp->lat = (1-(0.15 + index*(cp_rp->subframe_ID-1))) *
			  (0.25*a_lat + 0.75*d_lat) + 
			  (0.15 + index*(cp_rp->subframe_ID-1)) *
			  (0.25*b_lat + 0.75*e_lat); 
		rp->lon = (1-(0.15 + index*(cp_rp->subframe_ID-1))) *
			  (0.25*a_lon + 0.75*d_lon) + 
			  (0.15 + index*(cp_rp->subframe_ID-1)) *
			  (0.25*b_lon + 0.75*e_lon); 
	    } else if (cp_rp->subframe_ID > 4) {
		rp->lat = (1-(0.15 + index*(cp_rp->subframe_ID-5))) *
			  (0.5*a_lat + 0.5*d_lat) + 
			  (0.15 + index*(cp_rp->subframe_ID-5)) *
			  (0.5*b_lat + 0.5*e_lat); 
		rp->lon = (1-(0.15 + index*(cp_rp->subframe_ID-5))) *
			  (0.5*a_lon + 0.5*d_lon) + 
			  (0.15 + index*(cp_rp->subframe_ID-5)) *
			  (0.5*b_lon + 0.5*e_lon); 
	    }
	    if (rp->lon > 180) rp->lon -= 360;
	    printf("lat=%g,lon=%g\n",rp->lat,rp->lon); 

	} else {	
	    rp->lat = (float) ODLGetDouble(odl, CENTER_LAT, &status);
	    if ( status == FAIL ){
		asp_msg(LOG_DEBUG,"GetODLFrame: can't read %s", CENTER_LAT);
		return(FAIL);
	    }
	    rp->lon = (float) ODLGetDouble(odl, CENTER_LON, &status);
	    if ( status == FAIL ){
		asp_msg(LOG_DEBUG,"GetODLFrame: can't read %s", CENTER_LON);
		return(FAIL);
	    }
	}
	if (vbose) printf("lat=%g,lon=%g\n", rp->lat,rp->lon);
/*
	rp->lat = (float) ODLGetDouble(odl, CENTER_LAT, &status);
	if ( status == FAIL ){
		asp_msg(LOG_DEBUG,"GetODLFrame: can't read %s", CENTER_LAT);
		return(FAIL);
	}
	rp->lon = (float) ODLGetDouble(odl, CENTER_LON, &status);
	if ( status == FAIL ){
		asp_msg(LOG_DEBUG,"GetODLFrame: can't read %s", CENTER_LON);
		return(FAIL);
	}
*/
	return(PASS);
	
} /* end of GetODLFrame */

int ReadScanResults(seg_list,pre_list)
	TAPE_SEG_PTR seg_list;
	PREAM_FILE_PTR pre_list;
{
	int status;
	FILE *fd;
	char s[100], *str, tmp_file[100];
	int len;

  	strcpy(scan_result_filename," ");
   	strcpy(scan_result_version," ");
	use_rds_scan = 0;

        if ( cp_rp->proc_type == FRAME_RQST ){
	   strcpy(tmp_file,cp_rp->frame_file);
	   if (RcpCPToFile( cp_rp->frame_file, ASP_SCAN_FILE ) != PASS){
		asp_msg(LOG_ERR,asp_messages[RCP_ERR],cp_rp->frame_file);
		return(FAIL);
	   }

           /* extract the filename of scan_result from frame_request and
	      version # from scan_result file to put in the ceos_leader file */

	   if ((strstr(tmp_file,rev_name)) != NULL){
	        strcpy(scan_result_filename, strstr(tmp_file,rev_name));
                if (vbose) printf("scan_filename=%s\n", scan_result_filename);
	   } else printf("Error in getting scan_result_filename\n");

           if ( (fd = fopen(ASP_SCAN_FILE,"r")) == NULL ){
                asp_msg(LOG_DEBUG,"Error in opening file %s",ASP_SCAN_FILE);
                return(FAIL);
           }
           if (fgets(s,100,fd) == NULL){
                asp_msg(LOG_DEBUG,"Error in reading file %s",ASP_SCAN_FILE);
                return(FAIL);
           }
	
	   if ((strstr(s,"ASP")) == NULL) {
	    	if ((strstr(s,"RDS")) == NULL) {
		    printf("no version number in scan_result file\n");
		    strcpy(scan_result_version," ");
		}
		else {
		    printf("use RDS scan_result file\n");
		    str = strstr(s,"RDS");
	   	    len = strlen(str);
	   	    len -= 3;
	   	    strncpy(scan_result_version,str,len);
	   	    scan_result_version[len]='\0';
		    use_rds_scan = 1;	/* CV 7/24/97 */
		}
	   }
	   else {
		printf("use ASP scan_result file\n");
		str = strstr(s,"ASP");
	   	len = strlen(str);
	   	len -= 3;
	   	strncpy(scan_result_version,str,len);
	   	scan_result_version[len]='\0';
	   }
           if (vbose) printf("scan version=%s\n", scan_result_version);
	   fclose(fd);

           status = GetODLScanFile(ASP_SCAN_FILE, cp_rp->frame_ID,
                       seg_list, pre_list, Cur_Rqst);
	   if (status == FAIL) return(FAIL);
	   /* update job file */
	   write_job_file( Cur_Rqst, &sv1, &tcf );
	   return(status);
        } return(1);
} /* end of ReadScanResults */

int GetCPCalFile( data )
	float *data;
{
	if (RcpCPToFile( cp_rp->cal_file, ASP_PVS_FILE ) != PASS){
		asp_msg(LOG_ERR,asp_messages[RCP_ERR],cp_rp->cal_file);
		return(FAIL);
	}
	return( GetCalFile( ASP_PVS_FILE, data ) );
	
} /* end of GetCPCalFile */

int RcpCPToFile( infile, outfile )
	char *infile, *outfile;
{
	char cmd[200];
	char string[120],string1[120];

	strcpy(string,strtok(infile,":"));
	strcpy(string1,string);
	strcpy(string,strtok(NULL,"\0"));

	sprintf(cmd,"rcp %s-fddi:%s %s",string1,string,outfile);
	asp_msg(LOG_DEBUG,"rcp %s-fddi:%s %s\n",string1,string,outfile);
	return(system(cmd));

} /* end of RcpFileToCP */

int RcpFileToCP( infile, outfile )
	char *infile, *outfile;
{
	char cmd[132];
	char string[120],string1[120];

	strcpy(string,strtok(outfile,":"));
	strcpy(string1,string);
	strcpy(string,strtok(NULL,"\0"));
	sprintf(cmd,"rcp %s %s-fddi:%s",infile,string1,string);
	asp_msg(LOG_DEBUG,"rcp %s %s-fddi:%s\n",infile,string1,string);
	return(system(cmd));

} /* end of RcpFileToCP */

int RcpImageFile()
{
	if ( !GetResType() ){
	   if (RcpFileToCP(ASP_AVG_IMAGE_FILE, cp_rp->image_file) != PASS){
		asp_msg(LOG_ERR,asp_messages[RCP_ERR],cp_rp->image_file);
		return(FAIL);
	   }
	} else {
	   if (RcpFileToCP(ASP_FULL_IMAGE_FILE, cp_rp->image_file) != PASS){
		asp_msg(LOG_ERR,asp_messages[RCP_ERR],cp_rp->image_file);
		return(FAIL);
	   }
	}
	return(PASS);
} /* end of RcpImageFile */

int RcpLdrFile()
{
	if (RcpFileToCP(ASP_CEOS_LDR_FILE, cp_rp->ceos_ldr) != PASS){
	   asp_msg(LOG_ERR,asp_messages[RCP_ERR],cp_rp->ceos_ldr);
	   return(FAIL);
	}
	return(PASS);
} /* end of RcpLdrFile */


int WriteCleanupStatusToCP( msg, status, job_id )
	ODL msg;
	int status, job_id;
{
	ODL ack;

	if ((ack = GetAckMsg(msg, "SUBSYSTEM_STATUS")) == NULL) {
		asp_msg(LOG_DEBUG,"can't get SUBSYSTEM_STATUS\n");
		return(FAIL);
	}
	if (!ODLSetVal(ack,"SUBSYSTEM_STATUS.COMMON_HEADER.DESTINATION","CP")){
		asp_msg(LOG_DEBUG,"can't set COMMON.DESTINATION\n");
	   	return(FAIL);
	}
	if (!ODLSetVal(ack,"SUBSYSTEM_STATUS.COMMON_HEADER.SOURCE","ASP")){
		asp_msg(LOG_DEBUG,"can't set COMMON.SOURCE\n");
	   	return(FAIL);
	}
	if (!ODLSetVal(ack,"SUBSYSTEM_STATUS.BODY.SUB_TYPE","CLEANUP")){
		asp_msg(LOG_DEBUG,"can't set BODY.SUB_TYPE\n");
	   	return(FAIL);
	}
	if (!ODLSetVal(ack,"SUBSYSTEM_STATUS.BODY.JOB_ID",job_id)){
		asp_msg(LOG_DEBUG,"can't set BODY.JOB_ID\n");
	   	return(FAIL);
	}
	if ( status == PASS ){
/* ASP-CP 2.1, 11/21/96 */
	   if (!ODLSetVal(ack, "SUBSYSTEM_STATUS.BODY.STATUS", "COMPLETED")){
/*
	   if (!ODLSetVal(ack, "SUBSYSTEM_STATUS.BODY.STATUS", PASS)){
*/
		asp_msg(LOG_DEBUG,"can't set BODY.STATUS\n");
		return(FAIL);
	   }
	} else {
/* ASP-CP 2.1, 11/21/96 */
	   if (!ODLSetVal(ack,"SUBSYSTEM_STATUS.BODY.STATUS", "CANCEL/FAIL")){
/*
	   if (!ODLSetVal(ack,"SUBSYSTEM_STATUS.BODY.STATUS",FAIL)){
*/
		asp_msg(LOG_DEBUG,"can't set BODY.STATUS\n");
		return(FAIL);
	   }
	   if (!ODLSetVal(ack,"SUBSYSTEM_STATUS.BODY.COMMENT",asp_error_msg)){
		asp_msg(LOG_DEBUG,"can't set BODY.COMMENT\n");
	   	return(FAIL);
	   }
	} 
	WriteMsgToServer(ack);
	ODLFree(ack);
        Write_CP_Error_Status = FAIL;
        asp_error_msg[0] = '\0';
	return(PASS);

} /* end of WriteCleanupStatusToCP */

int WriteStatusToCP( msg, status )
	ODL msg;
	int status;
{
	ODL ack;

        if ((ack = GetAckMsg(msg, "SUBSYSTEM_COMPLETED")) == NULL) {
		asp_msg(LOG_DEBUG,"can't get SUBSYSTEM_COMPLETED\n");
		return(FAIL);
        }
	if ( status == PASS ){
/* ASP-CP 2.1, 11/21/96 */
	   if ( !ODLSetVal(ack, "SUBSYSTEM_COMPLETED.BODY.STATUS", "COMPLETED") ){
/*
	   if ( !ODLSetVal(ack, "SUBSYSTEM_COMPLETED.BODY.STATUS", PASS) ){
*/
		asp_msg(LOG_DEBUG,"can't set BODY.STATUS\n");
		return(FAIL);
	   }
	} else {
/* ASP-CP 2.1, 11/21/96 */
           if (status == -2) {
	     if ( !ODLSetVal(ack, "SUBSYSTEM_COMPLETED.BODY.STATUS", "WRONG_TAPE") ){
		asp_msg(LOG_DEBUG,"can't set BODY.STATUS\n");
		return(FAIL);
      	     }

           sprintf(asp_error_msg, "WRONG_TAPE: Request %s, Inlog %s", Cur_Rqst->tape_id, inlog);
           }
           else
	     if ( !ODLSetVal(ack, "SUBSYSTEM_COMPLETED.BODY.STATUS", "CANCEL/FAIL") ){
		asp_msg(LOG_DEBUG,"can't set BODY.STATUS\n");
		return(FAIL);
	     }
/*
	   if ( !ODLSetVal(ack, "SUBSYSTEM_COMPLETED.BODY.STATUS", status) ){
		asp_msg(LOG_DEBUG,"can't set BODY.STATUS\n");
		return(FAIL);
	   }
*/
	   if ( !ODLSetVal(ack, "SUBSYSTEM_COMPLETED.BODY.COMMENT",asp_error_msg) ){
		asp_msg(LOG_DEBUG,"can't set BODY.COMMENT\n");
	   	return(FAIL);
	   }
	} 
	WriteMsgToServer(ack);
	ODLFree(ack);
        Write_CP_Error_Status = FAIL;
        asp_error_msg[0] = '\0';
	return(PASS);

} /* end of WriteStatusToCP */

int WriteScanResults(seg_list,filename)
	TAPE_SEG_PTR seg_list;
	char *filename;
{
	FILE *fp;
	int i, seg_cnt, pp_cnt;
	time_t now;
	struct tm *gmt_now;
	TAPE_SEG_PTR sp;
	PP_BLOCK_PTR pp;
	ODL msg;
	char *tmp, *token, *GetODLTok();
	char str[8], *tmp1, *tmp2;
	int month, day, sec, mills;
	double dist2, dsec;
	int datatake_staddr, datatake_endaddr;
	GMT datatake_sttime, datatake_endtime;
	GMT center_gmt, gmt;
	void fprintf_array();
	int mode_mismatch;

	if ( (fp = fopen(filename,"w+")) == NULL ){
		asp_msg(LOG_DEBUG,"Error in opening file %s",filename);
		return(FAIL);
	}
	now = time(NULL);
	gmt_now = gmtime(&now);
	for (seg_cnt = 0, sp = seg_list; sp != NULL; 
			sp = sp->nxt_seg, seg_cnt++ ){
		if ( seg_cnt == 0 ){
			datatake_staddr = sp->blk_start;
			datatake_sttime = sp->start;
			datatake_endaddr = sp->blk_end;
			datatake_endtime = sp->end;
		} else { 
			datatake_endaddr = sp->blk_end;
			datatake_endtime = sp->end;
		}
	}

	fprintf(fp,"%s%s */\n\n",ASP_SCAN_RESULT_VERSION,ASP_VERSION);
	fprintf(fp,"OBJECT = SCAN_RESULTS_FILE\n");


/* Write COMMON_HEADER OBJECT */

	if (vbose) printf("writing COMMON_HEADER\n");

	fprintf(fp,"  OBJECT = COMMON_HEADER\n");
	fprintf(fp,"    TIME = %4.4d-%3.3dT%2.2d:%2.2d:%06.3f\n",
		1900+gmt_now->tm_year, gmt_now->tm_yday+1, gmt_now->tm_hour,
		gmt_now->tm_min, (float) gmt_now->tm_sec );
	fprintf(fp,"    MSG_TYPE = \"SCAN_RESULTS_FILE\"\n");
	fprintf(fp,"    DESTINATION = \"CP\"\n");
	fprintf(fp,"    SOURCE = \"ASP\"\n");
	fprintf(fp,"    NUMBER_OF_RECORDS = 1\n");
	fprintf(fp,"  END_OBJECT = COMMON_HEADER\n\n");

/* Write BODY OBJECT */

	if (vbose) printf("writing BODY\n");
	fprintf(fp,"  OBJECT = BODY\n");

	for ( i = 0; RqstToScan[i] != NULL; ++i )
		fprintf(fp,"    %s\n",GetODLTok(save_rqst,RqstToScan[i]));
	
	fprintf(fp,"    SEGMENT_COUNT = %d\n",seg_cnt);

        /* if no segment (no STD mode data), use the Cur_rqst information to
    	fill in the scan_result file */

	if (seg_cnt == 0) {
	   fprintf(fp,"    %s = %d\n",SEG_START_TAPE_BLOCK,Cur_Rqst->start_blk);
	   fprintf(fp,"    %s = %d\n",SEG_END_TAPE_BLOCK,Cur_Rqst->end_blk);
	   fprintf(fp,"    %s = %4.4d-%3.3dT%2.2d:%2.2d:%06.3f\n",
		SEG_START_TIME,
		Cur_Rqst->start.yr, Cur_Rqst->start.day, Cur_Rqst->start.hr,
		Cur_Rqst->start.min, Cur_Rqst->start.second );
	   fprintf(fp,"    %s = %4.4d-%3.3dT%2.2d:%2.2d:%06.3f\n",
		SEG_END_TIME,
		Cur_Rqst->end.yr, Cur_Rqst->end.day, Cur_Rqst->end.hr,
		Cur_Rqst->end.min, Cur_Rqst->end.second );
	} else {
	    fprintf(fp,"    %s = %d\n",SEG_START_TAPE_BLOCK,datatake_staddr);
	    fprintf(fp,"    %s = %d\n",SEG_END_TAPE_BLOCK,datatake_endaddr);
	    fprintf(fp,"    %s = %4.4d-%3.3dT%2.2d:%2.2d:%06.3f\n",
		SEG_START_TIME,
		datatake_sttime.yr, datatake_sttime.day, datatake_sttime.hr,
		datatake_sttime.min, datatake_sttime.second );
	    fprintf(fp,"    %s = %4.4d-%3.3dT%2.2d:%2.2d:%06.3f\n",
		SEG_END_TIME,
		datatake_endtime.yr, datatake_endtime.day, datatake_endtime.hr,
		datatake_endtime.min, datatake_endtime.second );
	}
	fprintf(fp,"    PRE_CAL1_POW = %f\n",0.0); /* TBD */
	fprintf(fp,"    PRE_CAL2_POW = %f\n",0.0); /* TBD */
	fprintf(fp,"    POST_CAL1_POW = %f\n",0.0); /* TBD */
	fprintf(fp,"    POST_CAL2_POW = %f\n",0.0); /* TBD */
        fprintf(fp,"    OBJECT = STATE_VECTOR_RECORD\n");

        msg = StrToODL(save_rqst, strlen(save_rqst)+1);
        tmp = (char *) ODLToStr(msg,STATE_VECTOR_METADATA);
        tmp1 = (char *) ODLToStr(msg,STATE_VECTOR_DATA);
        tmp2 = (char *) ODLToStr(msg,GHA_CORRECTION);
        ODLFree(msg);

        token = strtok(tmp,"\n");
        while ( token != NULL ){
           fprintf(fp,"       %s\n",token);
           token = strtok(NULL,"\n");
           if ( !strcmp(token,"END") ) token = NULL;
        }
        free(tmp);
	free(token);

        token = strtok(tmp1,"\n");
        while ( token != NULL ){
           fprintf(fp,"       %s\n",token);
           token = strtok(NULL,"\n");
           if ( !strcmp(token,"END") ) token = NULL;
        }
        free(tmp1);
	free(token);

        fprintf(fp,"    END_OBJECT = STATE_VECTOR_RECORD\n");

        token = strtok(tmp2,"\n");
        while ( token != NULL ){
           fprintf(fp,"    %s\n",token);
           token = strtok(NULL,"\n");
           if ( !strcmp(token,"END") ) token = NULL;
        }
        free(tmp2);
	free(token);

	for (seg_cnt = 1, sp = seg_list; sp != NULL; sp = sp->nxt_seg){ 

	/* Write SEGMENT OBJECT */
	   mode_mismatch = 0;
	   if (vbose) printf("writing SEGMENT \n");
	   if ( Cur_Rqst->take_id[0] == 'R' ) {
	       printf("mode in job_file =%s, scan show = %s\n", cp_rp->mode,
			R_mode[sp->aux.beam_seq]);
	       if (strcmp(R_mode[sp->aux.beam_seq],cp_rp->mode))
		  mode_mismatch = 1;
	   }
	   if (mode_mismatch) {
	       fprintf(fp,"/*    OBJECT = %s */\n",SEGMENT);
	       fprintf(fp,"/*      %s = \"%s\" */\n",SEG_MODE, 
				R_mode[sp->aux.beam_seq]);
               fprintf(fp,"/*      EXPERIMENTAL_BEAM = %d */\n", 
			experimental_beam);
	       fprintf(fp,"/*      SEGMENT_ID = %d */\n", seg_cnt++ );
	       fprintf(fp,"/*      %s = %d */\n",SEG_START_TAPE_BLOCK, 
			sp->blk_start);
	       fprintf(fp,"/*      %s = %d */\n",SEG_END_TAPE_BLOCK, 
			sp->blk_end);
	       fprintf(fp,"/*      %s = %4.4d-%3.3dT%2.2d:%2.2d:%06.3f */\n",
			SEG_START_TIME,
			sp->start.yr, sp->start.day, sp->start.hr,
			sp->start.min, sp->start.second );
	       fprintf(fp,"/*      %s = %4.4d-%3.3dT%2.2d:%2.2d:%06.3f */\n",
			SEG_END_TIME,
			sp->end.yr, sp->end.day, sp->end.hr,
			sp->end.min, sp->end.second );
	       fprintf(fp,"/*    END_OBJECT = SEGMENT */\n\n");
	       mode_mismatch = 0;
	   } else {
	      fprintf(fp,"    OBJECT = %s\n",SEGMENT);
	      if ( Cur_Rqst->take_id[0] == 'R' ) 
	        fprintf(fp,"      %s = \"%s\"\n",SEG_MODE, 
				R_mode[sp->aux.beam_seq]);
	      else
	        fprintf(fp,"      %s = \"%s\"\n",SEG_MODE, cp_rp->mode);
              fprintf(fp,"      EXPERIMENTAL_BEAM = %d\n", experimental_beam);
	      fprintf(fp,"      SEGMENT_ID = %d\n", seg_cnt++);
	      fprintf(fp,"      %s = %d\n",SEG_START_TAPE_BLOCK, sp->blk_start);
	      fprintf(fp,"      %s = %d\n",SEG_END_TAPE_BLOCK, sp->blk_end);
	      fprintf(fp,"      %s = %4.4d-%3.3dT%2.2d:%2.2d:%06.3f\n",
		SEG_START_TIME,
		sp->start.yr, sp->start.day, sp->start.hr,
		sp->start.min, sp->start.second );
	     fprintf(fp,"      %s = %4.4d-%3.3dT%2.2d:%2.2d:%06.3f\n",
		SEG_END_TIME,
		sp->end.yr, sp->end.day, sp->end.hr,
		sp->end.min, sp->end.second );
	     pp = sp->ppb_list;
	     for (pp_cnt = 0, pp = sp->ppb_list; pp != NULL; pp = pp->nxt_ppb) 
			++pp_cnt;

	     fprintf(fp,"      %s = %d\n",FRAME_COUNT,pp_cnt);

	     fprintf_array( fp, FORMAT_CHANGE, sp->win_count, sp->wfmt );
	     fprintf_array( fp, WINDOW_POSN_CHANGE, sp->win_count, sp->wdwp );
	     fprintf_array( fp, AGC_FORMAT_CHANGE, sp->agc_count, sp->afmt );
	     fprintf_array( fp, AGC_VALUE_CHANGE, sp->agc_count, sp->agc );

	     fprintf(fp,"      %s = %f\n",SCAN_PRF,
				sp->pre->prf); 
	     for (pp = sp->ppb_list; pp != NULL; pp = pp->nxt_ppb){
	     if (vbose) printf("writing FRAME \n");

	     /* Write FRAME OBJECT */

	     fprintf(fp,"      OBJECT = %s\n",IMAGE_FRAME);
	     fprintf(fp,"        %s = %d\n",FRAME_ID,pp->frame_id);
	     fprintf(fp,"        START_ADDRESS = %d\n",pp->blk_start);
	     fprintf(fp,"        END_ADDRESS = %d\n",pp->blk_end); 

	     fprintf(fp,"        START_TIME = %4.4d-%3.3dT%2.2d:%2.2d:%06.3f\n",
		pp->start_gmt.yr, pp->start_gmt.day, 
		pp->start_gmt.hr, pp->start_gmt.min, pp->start_gmt.second );

	     fprintf(fp,"        END_TIME = %4.4d-%3.3dT%2.2d:%2.2d:%06.3f\n",
		pp->end_gmt.yr, pp->end_gmt.day, 
		pp->end_gmt.hr, pp->end_gmt.min, pp->end_gmt.second );

/* These values can be obtained from scene_file structure for each region */

	     fprintf(fp,"        %s = %f\n",CENTER_LAT,pp->lat_rough); 
	     fprintf(fp,"        %s = %f\n",CENTER_LON,pp->lon_rough);
	     fprintf(fp, "        CENTER_TIME = %4.4d-%3.3dT%2.2d:%2.2d:%06.3f\n",
		pp->center_gmt.yr, pp->center_gmt.day, 
		pp->center_gmt.hr, pp->center_gmt.min, pp->center_gmt.second );
	     fprintf(fp,"        NEAR_START_LAT = %f\n",pp->lat_a); 
	     fprintf(fp,"        NEAR_START_LON = %f\n",pp->lon_a); 
	     fprintf(fp,"        NEAR_END_LAT = %f\n",pp->lat_d); 
	     fprintf(fp,"        NEAR_END_LON = %f\n",pp->lon_d);
	     fprintf(fp,"        FAR_START_LAT = %f\n",pp->lat_b); 
	     fprintf(fp,"        FAR_START_LON = %f\n",pp->lon_b); 
	     fprintf(fp,"        FAR_END_LAT = %f\n",pp->lat_e); 
	     fprintf(fp,"        FAR_END_LON = %f\n",pp->lon_e);
	     if ( pp->asc_dsc == asp_pass[ASC] )
	      fprintf(fp,"        ASC_DESC = \"%s\"\n",cp_pass[ASC]); 
	     else
	      fprintf(fp,"        ASC_DESC = \"%s\"\n",cp_pass[DSC]);

/* end scene_file */

	     fprintf(fp,"        START_FORMAT = %d\n",pp->fmt_start);
	     fprintf(fp,"        START_BIT_OFFSET = %d\n",pp->bit_off);
	     fprintf(fp,"        OBJECT = STATE_VECTOR_RECORD\n");

   	     msg = StrToODL(save_rqst, strlen(save_rqst)+1); 
	     tmp = (char *) ODLToStr(msg,STATE_VECTOR_METADATA);
	     ODLFree(msg);

	     token = strtok(tmp,"\n");
	     while ( token != NULL ){
	   	fprintf(fp,"          %s\n",token);
		token = strtok(NULL,"\n");
		if ( !strcmp(token,"END") ) token = NULL;
	     }
	     free(tmp);
	     free(token);

	     /*fprintf(fp,"          NUMBER_OF_VECTORS = %d\n",1);*/
	     fprintf(fp,"          OBJECT = STATE_VECTOR_DATA\n");
	     fprintf(fp,"            TIME = %4.4d-%3.3dT%2.2d:%2.2d:%06.3f\n",
		pp->sv.gmt.yr, pp->sv.gmt.day, pp->sv.gmt.hr,
		pp->sv.gmt.min, pp->sv.gmt.second );
	     fprintf(fp,"            X_POSITION = %f\n",pp->sv.pos.x);
	     fprintf(fp,"            Y_POSITION = %f\n",pp->sv.pos.y);
	     fprintf(fp,"            Z_POSITION = %f\n",pp->sv.pos.z);
	     fprintf(fp,"            X_VELOCITY = %f\n",pp->sv.vel.x);
	     fprintf(fp,"            Y_VELOCITY = %f\n",pp->sv.vel.y);
	     fprintf(fp,"            Z_VELOCITY = %f\n",pp->sv.vel.z);
	     fprintf(fp,"          END_OBJECT = STATE_VECTOR_DATA\n");

	     fprintf(fp,"        END_OBJECT = STATE_VECTOR_RECORD\n");
	   /* End FRAME OBJECT */

	     fprintf(fp,"      END_OBJECT = FRAME\n\n");

	     }	

	/* End SEGMENT OBJECT */

	     fprintf(fp,"    END_OBJECT = SEGMENT\n\n");
	  }
	}

/* End BODY OBJECT */
	fprintf(fp,"  END_OBJECT = BODY\n");

	fprintf(fp,"END_OBJECT = SCAN_RESULTS_FILE\n");
	fprintf(fp,"END\n");

	fclose(fp);
	return(PASS);
} /* end of WriteScanResults */

void fprintf_array( fp, KEY, count, array )
	FILE *fp;
	char *KEY;
	int  count;
	int  array[];
{
#define MAXITEM 8

	int i, cnt = count < MAXITEM ? count : MAXITEM;

	fprintf(fp,"      %s = (",KEY);
	for ( i = 0; i < cnt-1; i++ ) fprintf(fp,"%d, ",array[i]);
	if ( count > MAXITEM ){
		for ( i = MAXITEM-1; i < count-1; i++ ){
			if ( !(i % MAXITEM) ) fprintf(fp,"\n\t\t");
			fprintf(fp,"%d, ",array[i]);
		}
		if ( !(i % MAXITEM) ) fprintf(fp,"\n\t\t");
		fprintf(fp,"%d)\n",array[i]);
	} else if ( count >0 ) {
		  fprintf(fp,"%d)\n",array[i]);
	} else fprintf(fp,")\n");

} /* end of fprintf_array */

int WriteScanToCP(seg_list)
	TAPE_SEG_PTR seg_list;
{
	if ( WriteScanResults(seg_list,ASP_SCAN_FILE) == FAIL )
		return(FAIL);
	if (RcpFileToCP( ASP_SCAN_FILE, cp_rp->frame_file ) != PASS){
		asp_msg(LOG_ERR,asp_messages[RCP_ERR],cp_rp->frame_file);
		return(FAIL);
	}
	return(PASS);

} /* end of WriteScanToCP */

char *GetODLTok(char *source, char *keyword)
{
	char *tmp, *str;

	tmp = (char *) malloc(strlen(source)+1);
	strcpy(tmp,source);
	if ( strstr(tmp,keyword) != NULL ) {
	   str = strtok(strstr(tmp,keyword),"\n");
	} else {
	   printf("keyword: %s can not find \n",keyword);
	   str = NULL;
/*	   str = (char *) malloc(strlen(keyword)+16);
	   sprintf(str,"%s = %s",keyword,NULL); */
	}
	free(tmp);
	return(str);	
} /* end of GetODLTok */

int GetResType() 
{
	return( cp_rp->pix_spacing == cp_pix_spacing[LOW_RES] ? 
					LOW_RES : FULL_RES ); 
} /* end of GetResType */

int WritePMF(seg_list,filename)
	TAPE_SEG_PTR seg_list;
	char *filename;
{
	FILE *fp;
	time_t now;
	struct tm *gmt_now;
	char *GetODLTok();
	int type = RPR, factor = 1;
	register int i;
	char string[100], tmp_satellite[20], tmp_str[100];
	int len,len_prod_id;
	int is,fs;

	if ( (fp = fopen(filename,"w+")) == NULL ){
		asp_msg(LOG_DEBUG,"Error in opening file %s",filename);
		return(FAIL);
	}
	now = time(NULL);
	gmt_now = gmtime(&now);

	fprintf(fp,"OBJECT = %s\n",PMF_OBJ);

/* Write COMMON_HEADER OBJECT */

	if (vbose) printf("writing COMMON_HEADER\n");

	fprintf(fp,"  OBJECT = COMMON_HEADER\n");
        is = (int)(gmt_now->tm_sec);
        fs = (int)((gmt_now->tm_sec - (float)(is))*1000.0+0.5);
	fprintf(fp,"    TIME = %4.4d-%3.3dT%2.2d:%2.2d:%02.2d.%03.3d\n",
		1900+gmt_now->tm_year, gmt_now->tm_yday+1, gmt_now->tm_hour,
		gmt_now->tm_min, is, fs );
	fprintf(fp,"    MSG_TYPE = \"%s\"\n",PMF_OBJ);
	fprintf(fp,"    DESTINATION = \"CP\"\n");
	fprintf(fp,"    SOURCE = \"ASP\"\n");
	fprintf(fp,"    NUMBER_OF_RECORDS = 1\n");
	fprintf(fp,"  END_OBJECT = COMMON_HEADER\n\n");

/* Write CATALOG_METADATA OBJECT */

	if (vbose) printf("writing %s\n",PMF_META);
	fprintf(fp,"  OBJECT = %s\n",PMF_META);

	for ( i = 0; RqstToPMF[i] != NULL; ++i )
		fprintf(fp,"    %s\n",GetODLTok(save_rqst,RqstToPMF[i]));
	
	for ( i = 0; SegToPMF[i] != NULL; ++i )
		fprintf(fp,"    %s\n",GetODLTok(save_seg,SegToPMF[i]));
	
	for ( i = 0; FrameToPMF[i] != NULL; ++i )
		fprintf(fp,"    %s\n",GetODLTok(save_frame,FrameToPMF[i]));

	if (strncmp(Cur_Rqst->take_id,"E1",2) == 0)
	    strcpy(tmp_satellite,"ERS-1");
	else if (strncmp(Cur_Rqst->take_id,"E2",2) == 0)
	    strcpy(tmp_satellite,"ERS-2");
	else if (strncmp(Cur_Rqst->take_id,"J1",2) == 0)
	    strcpy(tmp_satellite,"JERS-1");
	else if (strncmp(Cur_Rqst->take_id,"R1",2) == 0)
	    strcpy(tmp_satellite,"RADARSAT-1");

	if ( ucs_product ) 
		sprintf(string,"%s",STD_BEAM_UNCOMP_PROD);
	else if (quicklook)
		sprintf(string,"%s",STD_BEAM_STD_PROD);
	else if ( !strcmp(Cur_Rqst->type, asp_product_type[CPX]) ) 
		sprintf(string,"%s",STD_BEAM_CPX_PROD);
	else if ( !strcmp(Cur_Rqst->type, asp_product_type[CSD]) ) 
		sprintf(string,"%s",STD_BEAM_CSD_PROD);
	else if ( !strcmp(Cur_Rqst->type, asp_product_type[RPR]) ) 
		sprintf(string,"%s",STD_BEAM_STD_PROD);

	fprintf(fp,"    DATASET = \"%s %s\"\n",tmp_satellite,string);
	free(string);
	fprintf(fp,"    FILE_VERSION = %0.3d\n",cp_rp->version); 

	len = strlen(cp_rp->pmf_file);
	len_prod_id = (strlen(cp_rp->product_id)) + 3;
	len = len - len_prod_id;
	strncpy(tmp_str,cp_rp->pmf_file,len);
	tmp_str[len] = 0;
	if (vbose) printf("string=%s,pmf_file=%s\n", tmp_str,cp_rp->pmf_file);	
	fprintf(fp,"    LOCAL_ARCHIVE_PATH = \"%s\"\n",tmp_str); /*constant*/
	if ( !strcmp(Cur_Rqst->type, asp_product_type[CPX]) ) 
		type = CPX;
	else if ( !strcmp(Cur_Rqst->type, asp_product_type[CSD]) )
		type = CSD;
	else if ( !GetResType() ) factor = 8;
		
	fprintf(fp,"    IMAGE_RECORD_COUNT = %d\n",
				(int)(asp_img_nrec[type]/factor+1) ); 
	fprintf(fp,"    IMAGE_FIRST_RECORD_LENGTH = %d\n",
				(int)(asp_img_reclen[type]/factor+192) );
	fprintf(fp,"    IMAGE_MAX_RECORD_LENGTH = %d\n",
				(int)(asp_img_reclen[type]/factor+192) );
	fprintf(fp,"    LEADER_MAX_RECORD_LENGTH = 5120\n"); 
	if ( !strcmp(Cur_Rqst->type, asp_product_type[CSD]) )
		fprintf(fp,"    LEADER_RECORD_COUNT = 7\n"); 
	else
		fprintf(fp,"    LEADER_RECORD_COUNT = 10\n");

	fprintf(fp,"    PRODUCT_CREATOR = \"ASP\"\n"); 
	fprintf(fp,"    CEOS_PRODUCT_TYPE = \"GRF\"\n"); /*funny*/

	fprintf(fp,"    PROC_VERSION = \"%s\"\n", ASP_VERSION); 
	fprintf(fp,"    STATE_VECTOR_PRECISION = \"%s\"\n",
				sv_precision[sv1.precision]); 
	if ( seg_list != NULL && seg_list->ppb_list != NULL )
	   fprintf(fp,"    SIGNAL_TO_NOISE_RATIO = %f\n",
					seg_list->ppb_list->snr);
	else fprintf(fp,"    SIGNAL_TO_NOISE_RATIO = \n");
	fprintf(fp,"    RADIOMETRIC_ACCURACY = 2.\n"); /*constant 2db*/
        is = (int)(gmt_now->tm_sec);
        fs = (int)((gmt_now->tm_sec - (float)(is))*1000.0+0.5);
	fprintf(fp,"    PRODUCT_CREATION_TIME = %4.4d-%3.3dT%2.2d:%2.2d:%02.2d.%03.3d\n",
		1900+gmt_now->tm_year, gmt_now->tm_yday+1, gmt_now->tm_hour,
		gmt_now->tm_min, is, fs);

/* End CATALOG_METADATA OBJECT */
	fprintf(fp,"  END_OBJECT = %s\n",PMF_META);
	fprintf(fp,"END_OBJECT = %s\n",PMF_OBJ);
	fprintf(fp,"END\n");

	fclose(fp);
	return(PASS);
} /* end of WritePMF */

int WritePMFToCP(seg_list)
	TAPE_SEG_PTR seg_list;
{
	if ( WritePMF(seg_list,ASP_PMF_FILE) == FAIL )
		return(FAIL);
	if (RcpFileToCP( ASP_PMF_FILE, cp_rp->pmf_file ) != PASS){
		asp_msg(LOG_ERR,asp_messages[RCP_ERR],cp_rp->pmf_file);
		return(FAIL);
	}
	return(PASS);

} /* end of WritePMFToCP */

int ReadCalFile( gain )
	float gain[900];
{
	char tmp_file[100];
	char id[10];

	strcpy(cal_params_filename," ");

	strcpy(tmp_file,cp_rp->cal_file);
	if (RcpCPToFile( cp_rp->cal_file, ASP_PVS_FILE ) != PASS){
		asp_msg(LOG_ERR,asp_messages[RCP_ERR],cp_rp->cal_file);
		return(FAIL);
	}
	sprintf(id,"%d",cp_rp->job_id);
	if (vbose) printf("tmp_file=%s,id=%s\n",tmp_file,id);
	if ((strstr(tmp_file,id)) != NULL){
	     strcpy(tmp_file,strstr(tmp_file,id));
	     strcpy(cal_params_filename,strtok(tmp_file,"/"));
	     strcpy(cal_params_filename,strtok(NULL,"\0"));
	     if (vbose) printf("cal_filename=%s\n",cal_params_filename);
	} else printf("Error in getting cal_params_filename\n");
 
	return( GetODLCalFile( ASP_PVS_FILE, gain ) );
	
} /* end of ReadCalFile */

int GetODLCalFile ( filename, data )
	char *filename;
	float data[900];
{
	ODL	odl,*gainvec,*Ant_Params, odl1;
	int i, j, nr = 1, nc = 900;
	char *string, src_antptn[20],file_name[100];
	int index, no_rec;
	double data_array[900], elev_incr, beamctr, first_elev;
	double aa_array[200],bb_array[200],a_array[200],b_array[200];
	double beam_inc,mid_value,x_prime;
	int x1, x2, total_no_gain, ctr_no;
	int start_pt, first_set, no_first_set;	
	int end_pt, last_set, no_last_set;	
	double start_angl, end_angl, scale_coef;
	int status = PASS;
	double total_gain;
	
	noise_fct = 0.0;
	linear_conv_fct = 0.0;
	offset_conv_fct = 0.0;
	islr_rng = pslr_rng = 0.0;
	range_res = azimuth_res = 0.0;
	range_ambig = azimuth_ambig = 0.0;
	ori_error = 0.0;
	dist_skew = 0.0;
	crosst_scale = alongt_scale = 0.0;
	crosst_locerr = alongt_locerr = 0.0;

	if ( (odl = GetODLFile( filename )) == NULL ){
	      asp_msg(LOG_DEBUG,"GetODLCalFile: can't read %s", filename);
	      ODLFree(odl); return(FAIL);
	}
	if ( Cur_Rqst->take_id[0] == 'R' ){
	        string = ODLGetString(odl, CAL_MODE, &status);
		if (status == FAIL) {
	            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",CAL_MODE);
		    free(string);
		    ODLFree(odl); return (FAIL);
		} else {
		    if (strcmp(string,cp_rp->mode)) {
	               asp_msg(LOG_DEBUG,"GetODLCalFile: incorrect beam_seq %s\n",string); 
		       ODLFree(odl);	
		       free(string);
	      	       return(FAIL);
		    }
		}	
	        if (vbose) 
		   printf("beam_seq=%s,cal_mode=%s\n", cp_rp->mode,string); 
		free(string);

	        string = ODLGetString(odl, SRC_ANTPTN, &status);
		if (status == FAIL) {
	           asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",SRC_ANTPTN);
		   free(string);
		   ODLFree(odl); return (FAIL);
		} else {
		   strncpy(src_antptn,string,3);
		   src_antptn[3] = '\0';
		   if (vbose) printf("src_antptn=%s\n",src_antptn); 
		   free(string);
		}
	}
        /* extract the cal_status & cal_comment from CPF */
        string = ODLGetString(odl, CAL_STATUS, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",CAL_STATUS);
            free(string);
            ODLFree(odl); return (FAIL);
        } else {
            strcpy(cal_status,string);
            free(string);
        }

        string = ODLGetString(odl, CAL_COMMENT, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",CAL_COMMENT);
            free(string);
            ODLFree(odl); return (FAIL);
        } else {
            strcpy(cal_comment,string);
            free(string);
        }
        if (vbose){
           printf("cal_status=%s\n",cal_status);
           printf("cal_comment=%s\n",cal_comment);
        }
        rel_radio_acc = ODLGetDouble(odl, REL_RADIO_ACC, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",REL_RADIO_ACC);
            ODLFree(odl); return (FAIL);
        }
        noise_fct = ODLGetDouble(odl, NOISE_FACT, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",NOISE_FACT);
            ODLFree(odl); return (FAIL);
        }
        linear_conv_fct = ODLGetDouble(odl, LINEAR_CONV_FACT, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",LINEAR_CONV_FACT);
            ODLFree(odl); return (FAIL);
        }
        offset_conv_fct = ODLGetDouble(odl, OFFSET_CONV_FACT, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",OFFSET_CONV_FACT);
            ODLFree(odl); return (FAIL);
        }
        islr_rng = ODLGetDouble(odl, CD_ISLR_RNG, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",CD_ISLR_RNG);
            ODLFree(odl); return (FAIL);
        }
        pslr_rng = ODLGetDouble(odl, CD_PSLR_RNG, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",CD_PSLR_RNG);
            ODLFree(odl); return (FAIL);
        }
        range_res = ODLGetDouble(odl, RNG_RES, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",RNG_RES);
            ODLFree(odl); return (FAIL);
        }
        azimuth_res = ODLGetDouble(odl, AZI_RES, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",AZI_RES);
            ODLFree(odl); return (FAIL);
        }
        range_ambig = ODLGetDouble(odl, RNG_AMBIG, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",RNG_AMBIG);
            ODLFree(odl); return (FAIL);
        }
        azimuth_ambig = ODLGetDouble(odl, AZI_AMBIG, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",AZI_AMBIG);
            ODLFree(odl); return (FAIL);
        }
        ori_error = ODLGetDouble(odl, ORI_ERR, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",ORI_ERR);
            ODLFree(odl); return (FAIL);
        }
        dist_skew = ODLGetDouble(odl, DIS_SKEW, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",DIS_SKEW);
            ODLFree(odl); return (FAIL);
        }
        crosst_scale = ODLGetDouble(odl, CRT_SCALE, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",CRT_SCALE);
            ODLFree(odl); return (FAIL);
        }
        alongt_scale = ODLGetDouble(odl, ALT_SCALE, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",ALT_SCALE);
            ODLFree(odl); return (FAIL);
        }
        crosst_locerr = ODLGetDouble(odl, CRT_LOCERR, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",CRT_LOCERR);
            ODLFree(odl); return (FAIL);
        }
        alongt_locerr = ODLGetDouble(odl, ALT_LOCERR, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",ALT_LOCERR);
            ODLFree(odl); return (FAIL);
        }

        beamctr = ODLGetDouble(odl, BEAMCTR, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",BEAMCTR);
            ODLFree(odl); return (FAIL);
        }
        elev_incr = ODLGetDouble(odl, ELEV_INCR, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",ELEV_INCR);
            ODLFree(odl); return (FAIL);
        }
        first_elev = ODLGetDouble(odl, FIRST_ELEV, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",FIRST_ELEV);
            ODLFree(odl); return (FAIL);
        }
        no_rec = ODLGetInt(odl, NO_REC, &status);
        if (status == FAIL) {
            asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",NO_REC);
            ODLFree(odl); return (FAIL);
        }

	if (ODLGetArrayDouble(odl,GAIN_VEC,data_array,&nr,&nc) == NULL){
		printf("GetODLCalFile: can't get %s\n",GAIN_VEC);
		ODLFree(odl); return(FAIL);
	}
	else ODLFree(odl);

	if (vbose) printf("nr=%d,nc=%d\n", nr,nc);
	for (i=0; i<nc; i++) {
	    data_array[i] = (pow(10.0,(-(data_array[i]/5.0))));	
/*	    printf("data_array[%d]=%g\n",i,data_array[i]); */
	}
	if ( Cur_Rqst->take_id[0] != 'R' ){
	    for (i=0; i<nc-1; i++) {
		data_array[i] = data_array[i+1];
	    }
	    normalize_data(data_array,nc-1);
	    for (i=0; i<nc-1; i++) {
		data[i] = (float)data_array[i];
	    }
	    ODLFree(odl);
	    return(nc-1);
	}
	else {
            sprintf(file_name,"%sant_params",PROC_PATH);
	    if ( (odl = GetODLFile( file_name )) == NULL ){
	       	asp_msg(LOG_DEBUG,"GetODLCalFile: can't read %s", file_name);
	       	ODLFree(odl); return(FAIL);
	    }
            if ((Ant_Params= (ODL*) Val(Lookup(odl,ANT_BODY))) == NULL){
               	asp_msg(LOG_DEBUG,"GetODLCalFile:can't read %s",file_name);
               	ODLFree(odl); return(FAIL);
            }
            for (i=0; Ant_Params[i]!=NULL; ++i) {
               	if (strcasecmp(Name(Ant_Params[i]), ANT_PARAMS) == 0) {
               	   string = ODLGetString(Ant_Params[i],ANT_MODE, &status);
               	   if (vbose) printf("string : %s\n",string);
               	   if (status == FAIL) {
               	      asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",ANT_MODE);
		      free(string);
           	      ODLFree(odl); return (FAIL);
                   } else {
                      if (!strcmp(cp_rp->mode,string)) {
                          free(string);
              	          start_angl = ODLGetDouble(Ant_Params[i],START_ANGLE,&status);
           	          if (status == FAIL) {
               	      	     asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",START_ANGLE);
               	      	     free(string);
               		     ODLFree(odl); return (FAIL);
		      	  }
           	          end_angl = ODLGetDouble(Ant_Params[i],END_ANGLE,&status);
           	          if (status == FAIL) {
               	             asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",END_ANGLE);
               	             free(string);
               	    	     ODLFree(odl); return (FAIL);
		          }
           	          total_gain = ODLGetDouble(Ant_Params[i],TOTAL_GAIN,&status);
           	          if (status == FAIL) {
               	      	     asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",TOTAL_GAIN);
               	      	     free(string);
               		     ODLFree(odl); return (FAIL);
		          }
           	          scale_coef = ODLGetDouble(Ant_Params[i],SCALE_COEF,&status);
           	          if (status == FAIL) {
               	      	     asp_msg(LOG_DEBUG,"GetODLCalFile:can't get %s\n",SCALE_COEF);
               	      	     free(string);
               		     ODLFree(odl); return (FAIL);
		          }
		          break;
		      } else
		          status = FAIL;
	           } 
                }
	    }  /* for loop */
	    if (status == FAIL) {
	        ODLFree(odl);
                asp_msg(LOG_DEBUG,"GetODLCalFile:Can not find info for %s\n", 			cp_rp->mode);
                free(string);
	   	return (FAIL);
	    }
	    ODLFree(odl);
	    if (!strcmp(src_antptn,"CSA")) { 
		beam_inc = 0.1/elev_incr;	
		ctr_no = (int)(no_rec/2)+1;
		mid_value = data_array[ctr_no-1];
		start_pt = nint((start_angl-first_elev)/elev_incr); 
		end_pt = nint((end_angl-first_elev)/elev_incr);
		no_first_set = ctr_no - start_pt;
		no_last_set = end_pt - ctr_no;
		first_set = nint(no_first_set/beam_inc);
		last_set = nint(no_last_set/beam_inc);
		total_no_gain = first_set * 2;

/* use the number in default file "ant_params" */
		total_no_gain = (int)(total_gain * 10);
		first_set = total_no_gain/2;

		if (vbose) {
		  printf("elev_inc=%g,beam_inc=%g,ctr_no=%d,mid_val=%d\n",
			elev_incr,beam_inc,ctr_no,mid_value);
		  printf("first:%d,%d,%d, ", start_pt,no_first_set,first_set);
		  printf("total_no_gain=%d\n",total_no_gain);
		}
		for (i=0; i<ctr_no-1; i++) {
		    a_array[i] = data_array[ctr_no+i];
		}
		for (i=0; i<ctr_no-1; i++) { 
		    b_array[i] = data_array[ctr_no-i-2];
		}
		for (i=0; i<first_set; i++) { /* set # of array to be same */ 
		    x_prime = 0.1/elev_incr * (i+1);
		    x1 = (int)((0.1/elev_incr) * (i+1));
		    x2 = x1+1;
		    aa_array[i] = a_array[x1] + ((a_array[x2] - a_array[x1])/
				  (x2 - x1)) * (x_prime - x1);
		    bb_array[i] = b_array[x1] + ((b_array[x2] - b_array[x1])/
				  (x2 - x1)) * (x_prime - x1);
		}
		for (i=0; i<total_no_gain; i++) {
		    if (i < first_set)
		       data_array[i] = bb_array[first_set-1-i];
		    else if (i == first_set)
		       data_array[i] = mid_value;
		    else
		       data_array[i] = aa_array[i-first_set-1];
		}				
/*	    	for (i=0; i<total_no_gain; i++) 
	        	printf("data[%d]=%10.7f\n",i,data_array[i]); 
*/
	    	normalize_data(data_array,total_no_gain);
 		for (i=0; i<total_no_gain; i++) {
			data[i] = (float)(data_array[i]*scale_coef);
/*	           	printf("data[%d]=%10.7f\n",i,data[i]); */
	    	}
	
	    	ODLFree(odl);
	    	return(total_no_gain);
	    } 
	    else {	 
	        for (i=0; i<nc-1; i++) 
			data_array[i] = data_array[i];
	    	normalize_data(data_array,nc-1);
	    	for (i=0; i<nc-1; i++) {
		    data[i] = (float)data_array[i]*scale_coef;
/*	            printf("data[%d]=%10.7f\n",i,data[i]); */
	    	}

	    	ODLFree(odl);
	    	return(nc-1);
	    }
	}
		
} /* end of GetODLCalFile */

normalize_data(data_array,cnt)
double *data_array;
int cnt;
{

        float large_no,normal_data[200];
        int i,j;
        FILE *fd3;

        if ((fd3 = fopen("gaintbl","w+")) == NULL){
                printf("Can not create the NORMALIZE file\n");
                return(FAIL);
        }
        j = 0;
        large_no = data_array[j];
        while (j < cnt){
          if (large_no < data_array[j]){
            large_no = data_array[j];
          }
          j += 1;
        }
        if (vbose) printf("large_no = %g\n", large_no);
	for (i=0; i<cnt; i++)
		normal_data[i]=data_array[i];	
        for (i=0; i<cnt; i++){
	   if ( Cur_Rqst->take_id[0] == 'R' )
                data_array[i] = data_array[i]/large_no;
	   else {
               if (normal_data[i] < large_no) {
                        data_array[i] = data_array[i]/large_no;
                } else
                        data_array[i]=large_no;
	   }
           fprintf(fd3,"%10.6f\n", data_array[i]);
        }
        fclose(fd3);
}

int GetODLFrameTable(beam_seq,filename,fix_ctr_tbl) 
char *beam_seq;
char *filename;
double *fix_ctr_tbl;
{
        ODL  odl, *FrameTable;
        int i, j, nr=1, nc=900;
        char *string, platform[10];
	int min_rev, max_rev;
	int status = PASS;

        if ( (odl = GetODLFile( filename )) == NULL ){
              asp_msg(LOG_DEBUG,"GetODLFrameTable: can't read %s", filename);
              ODLFree(odl); return(FAIL);
        }
	sprintf(platform,"%.2s", Cur_Rqst->take_id);
   	if ((FrameTable = (ODL*) Val(Lookup(odl,FRAME_TABLE_BODY))) == NULL){
              asp_msg(LOG_DEBUG,"GetODLFrameTable: can't read %s", filename);
      	      ODLFree(odl); return(FAIL);
   	}
	for (i=0; FrameTable[i]!=NULL; ++i) { 
         if (strcasecmp(Name(FrameTable[i]), FRAME_TABLE) == 0) {
	   string = ODLGetString(FrameTable[i],PLATFORM, &status);
	   if (vbose) printf("platform=%s,%s; ",platform,string);
           if (status == FAIL) {
               asp_msg(LOG_DEBUG,"GetODLFrameTable:can't get %s\n",PLATFORM);
               free(string);
               ODLFree(odl); return (FAIL);
	   }
	   else {
	      if (!strcmp(string,platform )) {
		  if ( !strcmp(string,"R1")){
		     free(string);
                     string = ODLGetString(FrameTable[i], TABLE_MODE, &status);
                     if (status == FAIL) {
                        asp_msg(LOG_DEBUG,"GetODLFrameTable:can't get %s\n",
				TABLE_MODE);
                        free(string);
                        ODLFree(odl); return (FAIL);
		     }
		     else {
                        if (strcmp(string,beam_seq)){
			   status = FAIL;
			   continue;
		        }
			free(string);	
		     }
		  }
                  min_rev = ODLGetInt(FrameTable[i], MIN_REV, &status);
                  if (status == FAIL) {
                     asp_msg(LOG_DEBUG,"GetODLFrameTable:can't get %s\n",
				MIN_REV);
                     ODLFree(odl); return (FAIL);
	          }
                  max_rev = ODLGetInt(FrameTable[i], MAX_REV, &status);
                  if (status == FAIL) {
                     asp_msg(LOG_DEBUG,"GetODLFrameTable:can't get %s\n",
				MAX_REV);
                     ODLFree(odl); return (FAIL);
	          }

		  if (vbose)
		     printf("min=%d,max=%d,rev=%d\n",min_rev,max_rev,sv1.rev);
	          if (sv1.rev >= min_rev && sv1.rev <= max_rev) {
                     if (ODLGetArrayDouble(FrameTable[i],CENTER_LATITUDES,fix_ctr_tbl,&nr,
			 &nc) ==  NULL){
                          printf("GetODLFrameTable: can't get %s\n",	
				CENTER_LATITUDES);
                          ODLFree(odl); return(FAIL);
		      }
		      status = PASS;
		      break;
                  } /* if sv1->rev */
		  else 
		     status = FAIL;
	       } /* if strcmp */
	       else
		   status = FAIL;
	    }
	  }
	} /* for loop */
        if (vbose) printf("nr=%d,nc=%d,status=%d\n", nr,nc,status);
/*        for (i=0; i<nc; i++) printf("%g; ",fix_ctr_tbl[i]); */
        ODLFree(odl);
        return(status);	
}

float GetPreCal1Pow (){ return(cal.pre_cal1_pow); }
float GetPreCal2Pow (){ return(cal.pre_cal2_pow); }
float GetPostCal1Pow (){ return(cal.post_cal1_pow); }
float GetPostCal2Pow (){ return(cal.post_cal2_pow); }

