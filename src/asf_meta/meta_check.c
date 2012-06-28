#include "asf.h"
#include "asf_meta.h"
#include "terrasar.h"
#include "radarsat2.h"
#include "airsar.h"
#include "xml_util.h"
#include <math.h>

int isCEOS(const char *dataFile, char **error)
{
  char **inBandName = NULL, **inMetaName = NULL, tmp[1024], *message = NULL;
  char baseName[512];
  int ret = TRUE, nBands, trailer;
  ceos_metadata_ext_t metadata_ext;
  ceos_data_ext_t     data_ext;

  if (strlen(dataFile) <= 0)
    return FALSE;

  char *dirName = (char *) MALLOC(sizeof(char)*1024);
  char *fileName = (char *) MALLOC(sizeof(char)*1024);
  split_dir_and_file(dataFile, dirName, fileName);
  metadata_ext = get_ceos_metadata_name(dataFile, &inMetaName, &trailer);
  data_ext = get_ceos_data_name(dataFile, baseName, &inBandName, &nBands);
  if (data_ext == NO_CEOS_DATA) {
    message = (char *) MALLOC(sizeof(char)*1024);
    strcpy(message, "");
    sprintf(tmp, "Data file%s of original data (%s) missing.\n",
	    (nBands>1) ? "s" : "", fileName);
    strcat(message, tmp);
    ret = FALSE;
  }    
  if (metadata_ext == NO_CEOS_METADATA) {
    if (!message) {
      message = (char *) MALLOC(sizeof(char)*1024);
      strcpy(message, "");
    }
    sprintf(tmp, "Metadata file of original data (%s) missing.\n", fileName);
    strcat(message, tmp);
    ret = FALSE;
  }

  *error = message;
  FREE(dirName);
  FREE(fileName);

  return ret;
}

static char *get_airsar(char *buf, char *str)
{
  char *p, *q;

  static char value[51];
  memset(value,0,51);

  q = (char *) CALLOC(51,sizeof(char));
  p = strstr(buf, str);
  if (p) {
    strncpy(q, p, 50);
    strcpy(value, q+strlen(str));
  }
  else
    strcpy(value, "");

  FREE(q);

  return value;
}

int isAIRSAR(char *dataFile)
{
  airsar_header *header = NULL;
  FILE *fp;
  int L_airsar = 0;
  int C_airsar = 0;
  int P_airsar = 0;
  char buf[4400], *value, *band_data, *s;
  double version;

  // Allocate memory and file handling
  value = (char *) MALLOC(sizeof(char)*25);
  header = (airsar_header *) CALLOC(1, sizeof(airsar_header));
  band_data = (char *) MALLOC(sizeof(char)*(strlen(dataFile)+32));
  strcpy(band_data, dataFile);

  // Try L-band
  s = strstr(band_data, "_");
  if (s) {
    *s = '\0';
    strcat(s, "_l.dat");
  }
  //if (fileExists(band_data)) {
    strcpy(buf, "");
    fp = fopen(band_data, "r");
    if (fp != NULL && fgets(buf, 4400, fp) == NULL)
      asfPrintError("Could not read general header\n");
    FCLOSE(fp);
    version = atof(get_airsar(buf, "JPL AIRCRAFT SAR PROCESSOR VERSION"));
    L_airsar = (version > 0.0) ? 1 : 0;
    /*}
  else
  L_airsar = FALSE;*/

  // Try C-band
  s = strstr(band_data, "_");
  if (s) {
    *s = '\0';
    strcat(s, "_c.dat");
  }
  //if (fileExists(band_data)) {
    strcpy(buf, "");
    fp = fopen(band_data, "r");
    if (fp != NULL && fgets(buf, 4400, fp) == NULL)
      asfPrintError("Could not read general header\n");
    FCLOSE(fp);
    version = atof(get_airsar(buf, "JPL AIRCRAFT SAR PROCESSOR VERSION"));
    C_airsar = (version > 0.0) ? 1 : 0;
    /*}
  else
  C_airsar = FALSE;*/

    // Try P-band
  s = strstr(band_data, "_");
  if (s) {
    *s = '\0';
    strcat(s, "_p.dat");
  }
  //if (fileExists(band_data)) {
    strcpy(buf, "");
    fp = fopen(band_data, "r");
    if (fp != NULL && fgets(buf, 4400, fp) == NULL)
      asfPrintError("Could not read general header\n");
    FCLOSE(fp);
    version = atof(get_airsar(buf, "JPL AIRCRAFT SAR PROCESSOR VERSION"));
    P_airsar = (version > 0.0) ? 1 : 0;
    /*}
  else
  P_airsar = FALSE;*/

  FREE(value);
  FREE(header);
  FREE(band_data);

  return (L_airsar || C_airsar || P_airsar);
}

int isTerrasar_ext(char *dataFile, int checkPolarimetry, char **error)
{
  int found = TRUE;
  // Let's first check for an .xml extension
  char *ext = findExt(dataFile);

  // If it has the correct extension, investigate it further
  // Might sound a little harsh but avoids some XML parser warning otherwise.
  if (ext && strcmp_case(ext, ".xml") == 0) {
    char *satellite = NULL;
    char tmp[256], *message = NULL, *path=NULL, *inDataName=NULL;
    char imageDataType[25];
    int ii, numberOfLayers;
    satellite = (char *) MALLOC(sizeof(char)*25);
    FILE *fp;
    fp = fopen(dataFile, "r");
    xmlDoc *doc = xmlReadFile(dataFile, NULL, 0);
    if (doc) {
      strcpy(satellite, xml_get_string_value(doc, 
        "level1Product.productInfo.missionInfo.mission"));

      // only care about TerraSAR-X data
      if (satellite && 
	  (strncmp_case(satellite, "TSX", 3) == 0 ||
	   strncmp_case(satellite, "TDX", 3) == 0)) {

	if (checkPolarimetry) {
	  strcpy(imageDataType, xml_get_string_value(doc, 
	     "level1Product.productInfo.imageDataInfo.imageDataType"));
	  if (strcmp_case(imageDataType, "COMPLEX") != 0) {
	    if (!message)
	      message = (char *) MALLOC(sizeof(char)*1024);
	    sprintf(tmp, "Wrong data type!\n"
		    "Polarimetric processing requires SLC data!\n");
	    strcat(message, tmp);
	    found = FALSE;
	  }
	}

	// path from the xml (metadata) file
	path = get_dirname(dataFile);
	inDataName = (char *) MALLOC(sizeof(char)*(strlen(path)+100));
	numberOfLayers = xml_get_int_value(doc, 
	  "level1Product.productInfo.imageDataInfo.numberOfLayers");

	for (ii=0; ii<numberOfLayers; ii++) {
	  if (strlen(path)>0) {
	    strcpy(inDataName, path);
	    if (inDataName[strlen(inDataName)-1] != '/')
	      strcat(inDataName, "/");
	  }
	  else
	    strcpy(inDataName, "");
	  
	  // check whether data file exists
	  strcat(inDataName, xml_get_string_value(doc, 
	    "level1Product.productComponents.imageData[%d].file.location.path",
	    ii));
	  strcat(inDataName, "/");
	  strcat(inDataName, xml_get_string_value(doc, 
	    "level1Product.productComponents.imageData[%d].file.location."
	    "filename", ii));
	  if (!fileExists(inDataName)) {
	    if (!message)
	      message = (char *) MALLOC(sizeof(char)*1024);
	    sprintf(tmp, "Data file (%s) does not exist!\n", inDataName);
	    strcat(message, tmp);
	    found = FALSE;
	  }	
	}
      }
      else
	found = FALSE;
    }
    fclose(fp);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    if (path)
      FREE(path);
    if (inDataName)
      FREE(inDataName);
  }

  return found;
}

int isTerrasar(char *dataFile, char **error)
{
  isTerrasar_ext(dataFile, FALSE, &error);
}

int isRadarsat2(char *dataFile, char **error)
{
  char dataType[25];
  int found = TRUE;
  // Let's first check for an .xml extension
  char *ext = findExt(dataFile);

  // If it has the correct extension, investigate it further
  // Might sound a little harsh but avoids some XML parser warning otherwise.
  if (ext && strcmp_case(ext, ".xml") == 0) {
    int ii, band_count = 0;
    char tmp[256], *path = NULL, *message = NULL, *inDataName = NULL;
    char polarizations[20];
    char *satellite = (char *) MALLOC(sizeof(char)*25);
    FILE *fp;
    fp = fopen(dataFile, "r");
    xmlDoc *doc = xmlReadFile(dataFile, NULL, 0);
    if (doc) {
      strcpy(satellite, 
	     xml_get_string_value(doc, "product.sourceAttributes.satellite"));
      
      // only care about Radarsat-2 data
      if (satellite &&
	  strcmp_case(satellite, "RADARSAT-2") == 0) {
	
	strcpy(dataType, xml_get_string_value(doc, 
          "product.imageAttributes.rasterAttributes.dataType"));	
	if (strcmp_case(dataType, "COMPLEX") != 0) {
	  if (!message)
	    message = (char *) MALLOC(sizeof(char)*1024);
	  sprintf(tmp, "Wrong data type!\n"
		  "Polarimetric processing requires SLC data!\n");
	  strcat(message, tmp);
	  found = FALSE;
	}
	
	// path from the xml (metadata) file
	path = get_dirname(dataFile);
	inDataName = (char *) MALLOC(sizeof(char)*(strlen(path)+100));
	strcpy(polarizations, xml_get_string_value(doc, 
	  "product.sourceAttributes.radarParameters.polarizations"));
	for (ii=0; ii<strlen(polarizations)-1; ii++)
	  if (polarizations[ii] == ' ')
	  polarizations[ii] = ',';
	if (strstr(polarizations, "HH"))
	  band_count++;
	if (strstr(polarizations, "VV"))
	  band_count++;
	if (strstr(polarizations, "HV"))
	  band_count++;
	if (strstr(polarizations, "VH"))
	  band_count++;
	
	for (ii=0; ii<band_count; ii++) {
	  if (strlen(path)>0) {
	    strcpy(inDataName, path);
	    if (inDataName[strlen(inDataName)-1] != '/')
	      strcat(inDataName, "/");
	  }
	  else
	    strcpy(inDataName, "");
	  strcat(inDataName, xml_get_string_value(doc, 
	    "product.imageAttributes.fullResolutionImageData[%d]", ii));
	  if (!fileExists(inDataName)) {
	    if (!message)
	      message = (char *) MALLOC(sizeof(char)*1024);
	    sprintf(tmp, "Data file (%s) does not exist!\n", inDataName);
	    strcat(message, tmp);
	    found = FALSE;
	  }
	}
      }
      else
	found = FALSE;
    }
    fclose(fp);
    xmlFreeDoc(doc);
    xmlCleanupParser();
  }    

  return found;
}

int isUAVSAR(char *dataFile, char **error)
{
  char dataType[25];
  int found = FALSE;
  // Let's first check for an .ann extension
  char *ext = findExt(dataFile);

  // If it has the correct extension, investigate it further
  if (ext && strcmp_case(ext, ".ann") == 0) {

    // We have the most generic annotation file ever, with no indication what
    // kind of data we are supposed to have. It just list everything that was
    // originally produced.
    // The only identifier for even UAVSAR, I could find, was the URL.
    char line[512];
    FILE *fp;
    fp = fopen(dataFile, "r");
    while (fgets(line, 512, fp)) {
      if (strstr(line, "uavsar.jpl.nasa.gov"))
	found = TRUE;
    }
    fclose(fp);
  }

  return found;
}

// Literally a copy of the lzFetch function. However, had to be done since
// the original functionality, written for Vexcel par files, indexed location
// blocks. This does not work with the location block in the metadata file.
int isspace(int c);

static char *getString(char *file,char *param, int *err)
{
  int structDepth, indexVec = 0;
  char structName[16][255];
  char line[255];
  FILE *in = fopen(file, "r");

  if (in == NULL) {
    if (err == NULL)
      asfPrintError("Couldn't open spec file '%s'!\n", file);
    else {
      *err = 3;
      return NULL;
    }
  }
    
  structDepth=0;
  strcpy(structName[structDepth],"");
  
  while (NULL != fgets(line, 255, in)) {
    if (strchr(line, '{') != NULL) {
      char newStruct[255];
      char index[255];
      sscanf(line, "%s", newStruct);
      if (0 == strcmp(newStruct, "state_vector")) {
	sprintf(index, "[%d]", indexVec);
	indexVec++;
      }
      else strcpy(index, "");
      strcpy(structName[structDepth+1], structName[structDepth]);
      if (structDepth != 0)
	strcat(structName[structDepth+1], ".");
      strcat(structName[structDepth+1], newStruct);
      strcat(structName[structDepth+1], index);
      structDepth++;
    } 
    else if (strchr(line,'}') != NULL)
      structDepth--;
    else {
      char thisParamName[255];
      char thisParam[255];
      sscanf(line, "%s", thisParamName);
      strcpy(thisParam, structName[structDepth]);
      if (structDepth != 0)
	strcat(thisParam, ".");
      strcat(thisParam, thisParamName);
      if (0 == strcmp(thisParam, param)) {
	char *ret = (char *) MALLOC(sizeof(char)*255);
	char *end = strchr(line, '#');
	end--;
	while (isspace(*end)) {
	  *end = '\0';
	  end--;
	}
	char *start = strchr(line, ':');
	start++;
	while (isspace(*start)) 
	  start++;
	strcpy(ret, start);
	if (err!=NULL)
	  *err=0; 
	FCLOSE(in);
	return ret;
      }
    }
  }
  if (err == NULL) {
    asfPrintWarning("Couldn't find field '%s' in spec file '%s'.\n",
		    param, file);
    FCLOSE(in);
    return NULL;
  }
  *err = 1;
  FCLOSE(in);
  return NULL;
}

static double getDouble(char *file, char *param, int *err)
{
  double ret;
  char *str;
  
  str = getString(file, param, err);
  if ((err != NULL) && (*err != 0)) 
    return 0.0;
  if (str == NULL && err != NULL) 
    *err = 2;
  if (str && 1 != sscanf(str, "%lf", &ret)) {
    if (err == NULL)
      asfPrintError("Couldn't convert '%s' to a double.\n",str);
    else
      *err=2;
  }
  FREE(str);
  return ret;
}

static int getInt(char *file, char *param, int *err)
{
    int ret;
    char *str = getString(file, param, err);
    if ((err != NULL) && (*err != 0)) 
      return 0;
    if (!str) 
      return 0;
    if (1 != sscanf(str, "%d", &ret)) {
      if (err==NULL) 
	asfPrintError("Couldn't convert '%s' to an integer.\n", str);
      else
	*err=3;
    }
    FREE(str);
    return ret;
}
