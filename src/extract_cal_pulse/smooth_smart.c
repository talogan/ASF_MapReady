void create_smoothed_dat(FILE *fpin_dat,char *filename)
{
  unsigned char buf[SAMPLES_PER_LINE];
  int **tbuf;
  int *tmp;
  int sum[SAMPLES_PER_LINE];
  unsigned char obuf[SAMPLES_PER_LINE];
  int curr_input_line, curr_output_line;
  int i,j;
  int cnt;
  FILE *fpout;
  
  curr_input_line = 0;
  curr_output_line = 0;
  for (j=0;j<SAMPLES_PER_LINE; j++) sum[j]=0;
  
  fpout = fopen(filename,"wb");
  tbuf = (int **) malloc (sizeof(int *)*KERNEL);
  for (i=0;i<KERNEL;i++) tbuf[i] = (int *) malloc(sizeof(int)*SAMPLES_PER_LINE);
  
  /* fill the kernel */
  for (i=0; i<KERNEL; i++) {
    cnt = fread(buf,1,SAMPLES_PER_LINE, fpin_dat);
    for (j=0;j<SAMPLES_PER_LINE; j++) {
      tbuf[curr_input_line][j] = (int) buf[j];
      sum[j] += tbuf[curr_input_line][j];
    }
    curr_input_line++;

    /* dump out first lines */
    if (i<KERNEL/2) { fwrite(buf,SAMPLES_PER_LINE,1,fpout); curr_output_line++; }
  }

  printf("entering kernel loop with cnt = %i\n",cnt);

  while (cnt == SAMPLES_PER_LINE) {
    
    /* write out the average of this kernel */
    for (j=0;j<SAMPLES_PER_LINE; j++) obuf[j] = (unsigned char) (sum[j] / KERNEL);
    fwrite(obuf,SAMPLES_PER_LINE,1,fpout); 
    curr_output_line++;
    
    /* read in the next line; remove old line; add in new line */
    cnt = fread(buf,1,SAMPLES_PER_LINE,fpin_dat);
    curr_input_line++;
    for (j=0;j<SAMPLES_PER_LINE; j++) sum[j] = sum[j] - tbuf[0][j] + (int)buf[j];
    
    /* swap pointers and refill the (new) last line in the kernel */
    tmp = tbuf[0];
    for (i = 0; i<KERNEL-1; i++) tbuf[i]=tbuf[i+1];
    tbuf[KERNEL-1] = tmp;
    for (j=0;j<SAMPLES_PER_LINE; j++) tbuf[KERNEL-1][j] = (int) buf[j];
    
  }
  
  /* dump out last lines */
  for (i=KERNEL/2+1; i<KERNEL; i++) {
    for (j=0;j<SAMPLES_PER_LINE; j++) obuf[j] = (unsigned char) tbuf[i][j];
    fwrite(obuf,SAMPLES_PER_LINE,1,fpout); 
    curr_output_line++;
  }
  
  fseek(fpin_dat,0,SEEK_SET);
  fclose(fpout);
  
}
