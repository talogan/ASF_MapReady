#include <ctype.h>
#include <asf.h>
#include <assert.h>

char *uc (const char *string)
{
  static char out[1024];

  int ii,n = strlen(string);
  if (n>1023) n=1023;

  for (ii=0; ii<n; ii++)
    out[ii]=toupper(string[ii]);
  out[ii]='\0';

  return out;
}

char *lc (const char *string)
{
  static char out[1024];

  int ii,n = strlen(string);
  if (n>1023) n=1023;

  for (ii=0; ii<n; ii++)
    out[ii]=tolower(string[ii]);
  out[ii]='\0';

  return out;
}

int strcmp_case(const char *s1, const char *s2)
{
    const char *p1 = s1;
    const char *p2 = s2;

    while (toupper(*p1)==toupper(*p2++))
        if (*p1++ == '\0') return 0;
    return (toupper(*p1) - toupper(*--p2));
}

int strncmp_case(const char *s1, const char *s2, size_t n)
{
    const char *p1 = s1;
    const char *p2 = s2;
    int i=0;

    while (toupper(*p1)==toupper(*p2++)) {
        if (*p1++ == '\0') return 0;
        if (++i == n) return 0;
    }
    return (toupper(*p1) - toupper(*--p2));
}

char *appendStr(const char *s1, const char *s2)
{
    char *ret = MALLOC((strlen(s1)+strlen(s2)+1)*sizeof(char));
    strcat(strcpy(ret, s1), s2);
    return ret;
}

// copies not more than len-1 characters from the string 'src'
// (including any terminating null characters), to 'dst'
// If strlen(src)<len, null characters are appended to dst
// to make the total len.
// If strlen(src)>=len, src is truncated to len-1 characters,
// and a null terminating character is appended to dst.
// Differs from strncpy in that:
//  - only copies len-1 characters, instead of len.
//  - dst is guaranteed null-terminated
char *strncpy_safe(char *dst, const char *src, size_t len)
{
    strncpy(dst, src, len-1);
    dst[len-1] = '\0';
    return dst;
}

char *trim_spaces(const char *s)
{
  char *dup = STRDUP(s);

  // strip end whitespace
  int n = strlen(dup)-1;
  while (isspace(dup[n]) && n>0)
    dup[n--] = '\0';

  // handle all-spaces case
  if (n==0)
    return STRDUP("");

  // now strip beginning whitespace, may assume not all spaces
  char *p = dup;
  while (isspace(*p))
    ++p;

  char *ret = STRDUP(p);
  FREE(dup);

  return ret;
}

void chomp(char *str)
{
  if (str[strlen(str)-1] == '\n')
    str[strlen(str)-1] = '\0';
}

// returns TRUE if "str" ends with the characters in "tail"
int endsWith(const char *str, const char *tail)
{
  if (strlen(tail) > strlen(str))
    return FALSE;
  else if (strlen(tail) == strlen(str))
    return strcmp_case(str, tail)==0;
  else
  {
    const char *p = str + strlen(str) - strlen(tail);
    return strcmp_case(p, tail)==0;
  }
}

// Reimplementation of strtok_r for Windows, mingw does not have it
char *STRTOK_R(char *str1, const char *str2, char **str3)
{
#ifdef win32
  char *ret;

  if (!str1)
    str1 = *str3;

  while (*str1 && strchr(str2, *str1))
    ++str1;
  if (*str1 == '\0')
    return NULL;

  ret = str1;

  while (*str1 && !strchr(str2, *str1))
    ++str1;

  if (*str1)
    *str1++ = '\0';

  *str3 = str1;

  return ret;
#else
  return strtok_r(str1, str2, str3);
#endif
}

// return the number of occurences of char 'c' in string 's'
int count_char(const char *s, char c)
{
  // allow passing in NULL for the string
  if (!s) return 0; 

  int i,n=0;
  for (i=0; i<strlen(s); ++i)
    if (s[i]==c) ++n;

  return n;
}

// a version of strstr() ignoring case.
char *strstr_case(const char *str, const char *key)
{
  const char *p1, *p2;

  while (*str != '\0') {
    p1 = str;
    p2 = key;
    while (tolower(*p2) == tolower(*p1) && *p1 != '\0') {
      p1++;
      p2++;
    }
    if (*p2 == '\0') {
      return (char *)str;
    }
    str++;
  }
  return NULL;
}

// replaces occurences in "str" of "searchStr" with "replaceStr"
char *strReplace(const char *str, const char *searchStr, const char *replaceStr)
{
  int searchStrLen = strlen(searchStr);
  int replaceStrLen = strlen(replaceStr);
  const char *curr=str, *next;
  int n=0;

  // first, just counting characters
  do {
    next = strstr(curr, searchStr);
    if (!next) {
      n += strlen(curr);
      break;
    }

    n += next-curr + replaceStrLen;
    curr = next + searchStrLen;
  }
  while (*curr != '\0');

  curr = str;
  char *ret = MALLOC(sizeof(char)*(n+2));
  strcpy(ret, "");

  // now, for reals
  do {
    next = strstr(curr, searchStr);
    if (!next) {
      strcat(ret, curr);
      break;
    }

    strncat(ret, curr, next-curr);
    strcat(ret, replaceStr);
    curr = next + searchStrLen;

    assert(strlen(ret) <= n);
  }
  while (*curr != '\0');

  return ret;
}

//---------------------------------------------------------------------------
// csv string parsing routines

static void strip_end_whitesp_inplace(char *s)
{
    char *p = s + strlen(s) - 1;
    while (isspace(*p) && p>s)
        *p-- = '\0';
}

static void consolidate_quotes(char *s)
{
    int i,j=0,l=strlen(s)+1;
    char *buf = MALLOC(sizeof(char)*l);
    for (i=0; i<l; ++i) {
        if (s[i]=='\"' && s[i+1]=='\"') {
            buf[j] = '\"';
            ++i;
        } else {
            buf[j] = s[i];
        }
        ++j;
    }
    strcpy(s,buf);
    FREE(buf);
}

char *quoted_string_parse(char *p, char *s, int max_len, int line_num)
{
  if (!p) {
    // reached end of the line prematurely
    if (line_num>0)
      printf("Line %d: Not enough data columns.\n", line_num);
    strcpy(s,"");
    return NULL;
  }

  // starting at p, eat characters, putting into s (allocated by caller),
  // until we reach the ending comma, or max_len.  For quoted strings,
  // we eat until we find the closing quote.

  // first, we eat whitespace
  while (isspace(*p))
    ++p;

  // now see if we have a quoted string or not
  if (*p == '\"') {
    // a quoted string... scan ahead to next quote instead of comma
    char *q = strchr(p+1, '\"');
    if (q) {
      // after the quote, the next non-whitespace character should be a comma
      // if it is another quote, then we need to keep going
      char *r = q+1;
      while (*r == '\"') {
        char *q1 = strchr(r+1, '\"');
        if (!q1)
          break;
        else {
          q = q1;
          r = q+1;
        }
      }

      // eat whitespace, until we get to what should be a comma
      while (isspace(*r))
        ++r;

      // now, eat the comma
      if (*r != ',') {
        if (line_num>0)
          printf("Line %d: Column entry has extra characters\n", line_num);
        // scan ahead to the next comma...
        r = strchr(r,',');        
      }
      
      // do not include the quotes in the returned value
      *q = '\0'; // temporary
      strncpy_safe(s, p+1, max_len);
      *q = '\"';

      // convert all instances of "" to just plain "
      consolidate_quotes(s);

      // now return pointer to just after the comma (if we found one)
      if (r)
        return r+1;
      else
        return NULL;
    }
    else {
      // no trailing quote found!
      // just put everything that's left into the output string
      strncpy_safe(s, p+1, max_len);
      // return NULL (to indicate this is the last column)
      return NULL;
    }
  }
  else {
    // non-quoted string
    char *q = strchr(p, ',');
    if (q) {
      *q = '\0'; // temporary
      strncpy_safe(s, p, max_len);
      strip_end_whitesp_inplace(s);
      *q = ',';

      // return pointer to just after the comma
      return q+1;
    }
    else {
      // no more commas, just return what is left in the string
      // and return NULL to indicate that we have everything now
      strncpy_safe(s, p, max_len);
      strip_end_whitesp_inplace(s);
      return NULL;
    } 
  }
}

// splits a quoted-element csv string into an array of elements
// for example, with:
//   str="hello","there ", you ," jerk!!"
//   sep=','
// you will get:
//   nelem=4
//   arr[0]=hello
//   arr[1]=there    <-- with a space at the end
//   arr[2]=you      <-- no spaces
//   arr[3]= jerk!!  <-- spaces preserved
//
// used like this:
//
//   char *str = ... read from file or something ...
//
//   int nelem;
//   char **arr;
//
//   split_into_array(str, ',', &nelem, &arr);
//
//   int i;
//   for (i = 0; i < nelem; ++i) {
//       // use arr[i]
//   }
//
//   free_char_array(&arr, nelem);
//
void split_into_array(char *str, char sep, int *nelem, char ***parr)
{
  char *p = str;
  int i,n=0;
  do {
    char val[512];
    p = quoted_string_parse(p, val, 512, -1);
    ++n;
  } while (p);

  char **arr = MALLOC(sizeof(char*)*n);
  p = str;
  for (i=0; i<n; ++i) {
    char val[512];
    p = quoted_string_parse(p, val, 512, -1);
    arr[i] = MALLOC(sizeof(char)*(strlen(val)+5));
    strcpy(arr[i], val);
  }

  *parr = arr;
  *nelem = n;
}

void free_char_array(char ***parr, int nelem)
{
  char **arr = *parr;
  int i;
  for (i=0; i<nelem; ++i)
    if (arr[i]) free(arr[i]);
  free(arr);
  *parr = NULL;
}

// splits a string into two pieces, stuff before the separater character
// and the stuff after it.  The separator character is not included in
// either string
void split2(const char *str_in, char sep, char **s1_out, char **s2_out)
{
  char *str = STRDUP(str_in);
  char *s1 = MALLOC(sizeof(char)*(strlen(str)+1));
  char *s2 = MALLOC(sizeof(char)*(strlen(str)+1));

  char *p = strchr(str, sep);

  if (p) {
    *p = '\0';
    strcpy(s1, str);
    *p = sep;
    strcpy(s2, p+1);
  } else {
    // no sep -- s2 is empty, s1 is a copy of str
    strcpy(s1, str);
    strcpy(s2, "");
  }

  // trim whitespace
  *s1_out = trim_spaces(s1);
  *s2_out = trim_spaces(s2);

  FREE(s1);
  FREE(s2);
  FREE(str);
}
