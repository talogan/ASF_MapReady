void create_smoothed_dat(FILE *fpin_dat,char *filename)
{
  unsigned char buf[SAMPLES_PER_LINE];
  int **tbuf;
  int *tmp;
  int sum[SAMPLES_PER_LINE];
  unsigned char obuf[SAMPLES_PER_LINE];
  int curr_input_line, curr_output_line, curr_line;
  int i,j;
  int cnt;
  FILE *fpout;
  
  curr_input_line = 0;
  curr_output_line = 0;
  curr_line = 0;
  
  fpout = fopen(filename,"wb");

  /* write out first lines */
  for (i=0;i<KERNEL/2; i++) {
    cnt = fread(buf,1,SAMPLES_PER_LINE, fpin_dat);
    fwrite(buf,1,SAMPLES_PER_LINE, fpout);
    curr_line++;
  }
  
  while (cnt == SAMPLES_PER_LINE) {
    for (j=0;j<SAMPLES_PER_LINE; j++) sum[j] = 0;
    for (line = curr_line-KERNEL/2; line < curr_line+KERNEL/2+1; line++)
      {
        long int where = (long int) (line*SAMPLES_PER_LINE);
        fseek(fpin_dat,where,SEEK_SET);
	cnt = fread(buf,1,SAMPLES_PER_LINE, fpin_dat);
        for (j=0;j<SAMPLES_PER_LINE; j++) sum[j] += (int) buf[j];
      }
    for (j=0;j<SAMPLES_PER_LINE; j++) obuf[j] = (unsigned int) (sum[j]/KERNEL);
    fwrite(obuf,1,SAMPLES_PER_LINE, fpout);
    curr_line++;
  }
  
  /* dump out last lines */
  for (i=KERNEL/2+1; i<KERNEL; i++) {
    long int where = (long int) (curr_line*SAMPLES_PER_LINE);
    fseek(fpin_dat,where,SEEK_SET);
    cnt = fread(buf,1,SAMPLES_PER_LINE, fpin_dat);
    fwrite(buf,1,SAMPLES_PER_LINE, fpout);
    curr_line++;
  }
  
  fseek(fpin_dat,0,SEEK_SET);
  fclose(fpout);
  
}
