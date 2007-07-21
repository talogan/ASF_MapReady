#include "asf_view.h"

typedef struct {
    FILE *fp;    // data file pointer
    int is_rgb;  // are we doing rgb compositing
    int band_gs; // which band we are using (when viewing as greyscale)
    int band_r;  // which band we are using for red (when rgb compositing)
    int band_g;  // which band we are using for green (when rgb compositing)
    int band_b;  // which band we are using for blue (when rgb compositing)
} ReadAsfClientInfo;

int try_asf(const char *filename)
{
    char *ext = findExt(filename);

    if (ext && strlen(ext) > 0) {
        return strcmp_case(ext, ".img") == 0 ||
               strcmp_case(ext, ".meta") == 0;
    } else {
        return try_ext(filename, ".img");
    }
}

int handle_asf_file(const char *filename, char *meta_name, char *data_name,
                    char **err)
{
    char *ext = findExt(filename);
    int has_ext = ext && strlen(ext) > 0;
    int has_asf_ext = has_ext &&
        (strcmp_case(ext,".img")==0 || strcmp_case(ext,".meta")==0);

    // either they gave us an ASF Internal extension, or we try adding it
    // to a user-provided basename, and that file does exist
    if (has_asf_ext || try_ext(filename, ".img"))
    {
        char *m = appendExt(filename, ".meta");
        char *d = appendExt(filename, ".img");

        strcpy(meta_name, m);
        strcpy(data_name, d);

        int ret;
        if (!fileExists(meta_name) || !fileExists(data_name)) {
            int l = sizeof(char)*strlen(filename)*2+255;
            *err = MALLOC(l);
            snprintf(*err, l,
                "Error opening ASF Internal Format file.\n"
                "  Metadata file: %s - %s\n"
                "      Data file: %s - %s\n",
                m, fileExists(m) ? "Found" : "NOT FOUND",
                d, fileExists(d) ? "Found" : "NOT FOUND");

            ret = FALSE;
        }
        else
            ret = TRUE;

        free(m);
        free(d);

        return ret;
    }
    else {
        // in theory this shouldn't happen, if try_ext is working
        assert(!try_ext(filename, ".img"));
        int l = sizeof(char)*strlen(filename)*2+255;
        *err = MALLOC(l);
        snprintf(*err, l,
            "Failed to open %s as an ASF Internal Format File.\n", filename);
        return FALSE;
    }

    // not reached
    assert(FALSE);
    return FALSE;
}

meta_parameters *read_asf_meta(const char *meta_name)
{
    return meta_read(meta_name);
}

int read_asf_client(int row_start, int n_rows_to_get,
                    void *dest_void, void *read_client_info,
                    meta_parameters *meta)
{
    ReadAsfClientInfo *info = (ReadAsfClientInfo*) read_client_info;
    int nl = meta->general->line_count;
    int ns = meta->general->sample_count;

    if (meta->general->data_type == BYTE) {
        unsigned char *dest = (unsigned char*)dest_void;
        if (data_ci->data_type==GREYSCALE_BYTE) {
            // reading byte data directly into the byte cache
            FSEEK64(info->fp, ns*(row_start + nl*info->band_gs), SEEK_SET);
            FREAD(dest, sizeof(unsigned char), n_rows_to_get*ns, info->fp);
        }
        else {
            // Here we have to read three separate strips of the
            // file to compose into the byte cache (which has interleaved
            // rgb values -- red vals coming from the first strip we read,
            // green from the second, and blues from the third.
            // So, we need a temporary buffer to place the values, so they
            // can be interleaved (i.e., we can't read directly into the
            // cache's memory)
            unsigned char *buf = MALLOC(sizeof(unsigned char)*ns);
            
            // first set the cache's memory to all zeros, this way any
            // rgb channels we don't populate will end up black
            memset(dest, 0, n_rows_to_get*ns*3);

            // red
            if (info->band_r >= 0) {
                int i,j,off = ns*(row_start + nl*info->band_r);
                for (i=0; i<n_rows_to_get; ++i) {
                    int k = 3*ns*i;
                    FSEEK64(info->fp, off + i*ns, SEEK_SET);
                    FREAD(buf, sizeof(unsigned char), ns, info->fp);
                    for (j=0; j<ns; ++j, k += 3)
                        dest[k] = buf[j];
                }
            }

            // green
            if (info->band_g >= 0) {
                int i,j,off = ns*(row_start + nl*info->band_g);
                for (i=0; i<n_rows_to_get; ++i) {
                    int k = 3*ns*i+1;
                    FSEEK64(info->fp, off + i*ns, SEEK_SET);
                    FREAD(buf, sizeof(unsigned char), ns, info->fp);
                    for (j=0; j<ns; ++j, k += 3)
                        dest[k] = buf[j];
                }
            }

            // blue
            if (info->band_b >= 0) {
                int i,j,off = ns*(row_start + nl*info->band_b);
                for (i=0; i<n_rows_to_get; ++i) {
                    int k = 3*ns*i+2;
                    FSEEK64(info->fp, off + i*ns, SEEK_SET);
                    FREAD(buf, sizeof(unsigned char), ns, info->fp);
                    for (j=0; j<ns; ++j, k += 3)
                        dest[k] = buf[j];
                }
            }

            free(buf);
        }
    } else {
        // this is the normal case, just reading in a strip of lines
        // from the file directly into the floating point cache
        assert(data_ci->data_type==GREYSCALE_FLOAT);
        float *dest = (float*)dest_void;
        get_float_lines(info->fp, meta, row_start + nl*info->band_gs,
                        n_rows_to_get, dest);
    }

    return TRUE;
}

int get_asf_thumbnail_data(int thumb_size_x, int thumb_size_y,
                           meta_parameters *meta, void *read_client_info,
                           void *dest_void)
{
    ReadAsfClientInfo *info = (ReadAsfClientInfo*) read_client_info;

    int sf = meta->general->line_count / thumb_size_y;
    assert(sf==meta->general->sample_count / thumb_size_x);
    int i,j;

    int nl = meta->general->line_count;
    //int ns = meta->general->sample_count;

    float *buf = MALLOC(sizeof(float)*meta->general->sample_count);

    if (meta->general->data_type == BYTE) {
        // BYTE case -- data file contains bytes.
        unsigned char *dest = (unsigned char*)dest_void;
        if (data_ci->data_type == GREYSCALE_BYTE) {
            // data file contains byte data, and we are just pulling out
            // one band to display.
            int off = nl*info->band_gs;
            for (i=0; i<thumb_size_y; ++i) {
                get_float_line(info->fp, meta, i*sf + off, buf);
                for (j=0; j<thumb_size_x; ++j)
                    dest[i*thumb_size_x+j] = (unsigned char)(buf[j*sf]);
                asfPercentMeter((float)i/(thumb_size_y-1));
            }
        }
        else {
            // rgb case -- we have to read up to 3 bands from the file,
            // and put them together
            assert(data_ci->data_type == RGB_BYTE);
            int off, k;

            // first set dest buffer to all zeros.  This way, any bands that
            // aren't set up will just come out black
            memset(dest, 0, 3*thumb_size_x*thumb_size_y);

            // "tot" is to help with the PercentMeter -- the total number
            // of lines that we will need to read.  "l" is the counter
            int tot = (info->band_r>=0) + (info->band_g>=0) + (info->band_b>=0);
            tot *= thumb_size_y;
            int l=0;

            printf("%d %d\n", tot, thumb_size_y);
            // to do each of the bands, we read the data into a float array,
            // then cast (back) to byte into the interleaved "dest" array
            // (interleaved in the sense that we only populate every 3rd item
            // each time through)

            // red band
            if (info->band_r >= 0) {
                off = nl*info->band_r;
                for (i=0; i<thumb_size_y; ++i) {
                    k=3*i*thumb_size_x; // starting point in dest array
                    get_float_line(info->fp, meta, i*sf + off, buf);
                    for (j=0; j<thumb_size_x; ++j, k += 3)
                        dest[k] = (unsigned char)(buf[j*sf]);
                    asfPercentMeter((float)(l++)/(tot-1));
                }
            }

            // green band
            if (info->band_g >= 0) {
                off = nl*info->band_g;
                for (i=0; i<thumb_size_y; ++i) {
                    k=3*i*thumb_size_x+1;
                    get_float_line(info->fp, meta, i*sf + off, buf);
                    for (j=0; j<thumb_size_x; ++j, k += 3)
                        dest[k] = (unsigned char)(buf[j*sf]);
                    asfPercentMeter((float)(l++)/(tot-1));
                }
            }

            // blue band
            if (info->band_b >= 0) {
                off = nl*info->band_b;
                for (i=0; i<thumb_size_y; ++i) {
                    k=3*i*thumb_size_x+2;
                    get_float_line(info->fp, meta, i*sf + off, buf);
                    for (j=0; j<thumb_size_x; ++j, k += 3)
                        dest[k] = (unsigned char)(buf[j*sf]);
                    asfPercentMeter((float)(l++)/(tot-1));
                }
            }

            //assert(l==tot);
            if (l!=tot) printf("%d %d\n", l, tot-1);
        }
    } else {
        // this is the normal case -- regular old floating point data,
        // we just read with get_float_line and populate directly into
        // a floating point array
        assert(data_ci->data_type==GREYSCALE_FLOAT);
        int off = nl*info->band_gs;
        float *dest = (float*)dest_void;
        for (i=0; i<thumb_size_y; ++i) {
            get_float_line(info->fp, meta, i*sf + off, buf);
            for (j=0; j<thumb_size_x; ++j)
                dest[i*thumb_size_x+j] = buf[j*sf];
            asfPercentMeter((float)i/(thumb_size_y-1));
        }
    }

    free(buf);
    return TRUE;
}

void free_asf_client_info(void *read_client_info)
{
    ReadAsfClientInfo *info = (ReadAsfClientInfo*) read_client_info;
    if (info->fp) fclose(info->fp);
    free(info);
}

int open_asf_data(const char *filename, const char *band,
                  meta_parameters *meta, ClientInterface *client)
{
    ReadAsfClientInfo *info = MALLOC(sizeof(ReadAsfClientInfo));

    info->is_rgb = FALSE;
    info->band_gs = info->band_r = info->band_g = info->band_b = 0;

    if (band) {
        char *r, *b, *g;
        if (split3(band, &r, &g, &b, ',')) {
            // Looks like we were given 3 bands -- so, we are doing rgb
            info->band_r = get_band_number(meta->general->bands,
                    meta->general->band_count, r);
            if (info->band_r < 0)
                asfPrintWarning("Red band '%s' not found.\n", r);
            else
                asfPrintStatus("Red band is band #%d: %s\n",
                    info->band_r+1, r);

            info->band_g = get_band_number(meta->general->bands,
                    meta->general->band_count, g);
            if (info->band_g < 0)
                asfPrintWarning("Green band '%s' not found.\n", g);
            else
                asfPrintStatus("Green band is band #%d: %s\n",
                    info->band_g+1, g);

            info->band_b = get_band_number(meta->general->bands,
                    meta->general->band_count, b);
            if (info->band_b < 0)
                asfPrintWarning("Blue band '%s' not found.\n", b);
            else
                asfPrintStatus("Blue band is band #%d: %s\n",
                    info->band_b+1, b);

            if (info->band_r < 0 && info->band_g < 0 && info->band_b < 0) {
                // none of the bands were found
                return FALSE;
            }

            info->is_rgb = TRUE;
            FREE(r); FREE(g); FREE(b);

            set_bands_rgb(info->band_r, info->band_g, info->band_b);
        } else {
            // Single band name given
            info->band_gs = get_band_number(meta->general->bands,
                    meta->general->band_count, (char*)band);
            if (info->band_gs < 0) {
                asfPrintWarning("Band '%s' not found.\n", band);
                return FALSE;
            } else {
                asfPrintStatus("Reading band #%d: %s\n",
                    info->band_gs+1, band);
            }

            set_bands_greyscale(info->band_gs);
        }
    }


    info->fp = fopen(filename, "rb");
    if (!info->fp) {
        asfPrintWarning("Failed to open ASF Internal file %s: %s\n",
            filename, strerror(errno));
        return FALSE;
    }

    client->read_client_info = info;
    client->read_fn = read_asf_client;
    client->thumb_fn = get_asf_thumbnail_data;
    client->free_fn = free_asf_client_info;

    if (meta->general->data_type == BYTE)
        client->data_type = info->is_rgb ? RGB_BYTE : GREYSCALE_BYTE;
    else
        client->data_type = GREYSCALE_FLOAT;

    //client->require_full_load = TRUE;
    return TRUE;
}