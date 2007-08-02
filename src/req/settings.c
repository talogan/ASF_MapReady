#include "req.h"
#include <ctype.h>

Settings *settings_new()
{
    Settings *s = MALLOC(sizeof(Settings));
    s->csv_dir = NULL;
    s->output_dir = NULL;
    s->req_num = 1;
    s->req_id = 1;
    s->station_code = NULL;
    return s;
}

void settings_free(Settings *s)
{
    if (s->csv_dir) FREE(s->csv_dir);
    if (s->output_dir) FREE(s->output_dir);
    if (s->station_code) FREE(s->station_code);
    FREE(s);
}

static int matches(const char *buf, const char *key)
{
    return strncmp_case(buf, key, strlen(key)) == 0;
}

static void read_string_param(const char *buf, const char *key, char **value)
{
    if (matches(buf, key)) {
        char *p = strchr(buf, '=');

        // skip past = sign, eat whitespace
        ++p;
        while (isspace(*p))
            ++p;

        *value = MALLOC(sizeof(char)*(strlen(p)+1));
        strcpy(*value, p);

        // eat trailing whitespace, too
        p = *value + strlen(*value) - 1;
        while (isspace(*p))
            *p-- = '\0';
    }
}

static void read_int_param(const char *buf, const char *key, int *value)
{
    if (matches(buf, key)) {
        char *p = strchr(buf, '=') + 1;
        *value = atoi(p);
    }
}

Settings *settings_load()
{
    Settings *s = settings_new();

    char *sav_file = find_in_share("req_settings.txt");
    if (fileExists(sav_file)) {
        FILE *fp = FOPEN(sav_file, "r");
        if (fp) {
            char buf[1024];
            while (fgets(buf, 1024, fp) != NULL) {
                read_string_param(buf, "csv directory", &s->csv_dir);
                read_string_param(buf, "output directory", &s->output_dir);
                read_int_param(buf, "next request number", &s->req_num);
                read_int_param(buf, "next request id", &s->req_id);
                read_string_param(buf, "station code", &s->station_code);
            }
            FCLOSE(fp);
        } else {
            printf("Failed to open settings file: %s\n", sav_file);
        }
    } else {
        printf("No settings file found.\n");
    }

    FREE(sav_file);

    if (s->req_num < 1) s->req_num = 1;
    if (s->req_id < 1) s->req_id = 1;

    if (!s->station_code)
        s->station_code = STRDUP("AADN");

    return s;
}

static void apply_settings_to_gui(Settings *s)
{
    put_string_to_entry("csv_directory_entry", s->csv_dir);
    put_string_to_entry("output_directory_entry", s->output_dir);
    put_int_to_entry("next_request_number_entry", s->req_num);
    put_int_to_entry("next_request_id_entry", s->req_id);
}

void apply_saved_settings()
{
    Settings *s = settings_load();
    apply_settings_to_gui(s);
    settings_free(s);
}

static void add_to_text(char **txt, int *cur_len, const char *format, ...)
{
    char buf[1024];

    va_list ap;
    va_start(ap, format);
    int len = vsnprintf(buf, sizeof(buf), format, ap);
    va_end(ap);

    if (len > 1020)
        printf("Warning! lengthy line may have been truncated.\n");

    if (strlen(*txt) + len >= *cur_len) {
        *cur_len += 1024;
        char *new_txt = MALLOC(sizeof(char)*(*cur_len));
        strcpy(new_txt, *txt);
        strcat(new_txt, buf);
        free(*txt);
        *txt = new_txt;
    } else {
        strcat(*txt, buf);
    }
}

void settings_save(Settings *s)
{
    char *sav_file = find_in_share("req_settings.txt");
    if (sav_file && fileExists(sav_file)) {
        // update an existing file
        printf("Found settings file: %s\n", sav_file);
        char *new_sav_txt = MALLOC(sizeof(char)*1024);
        strcpy(new_sav_txt, "");
        int len=1024;
        int wrote_csv=FALSE, wrote_output=FALSE, wrote_req_num=FALSE,
            wrote_req_id=FALSE, wrote_station_code=FALSE;
        FILE *fp = FOPEN(sav_file, "r");
        if (!fp) {
            message_box("Error opening output file!\n");
            return;
        }

        char buf[1024];
        while (fgets(buf, 1024, fp) != NULL) {
            if (s->csv_dir && matches(buf, "csv directory")) {
                add_to_text(&new_sav_txt, &len,
                    "csv directory = %s\r\n", s->csv_dir);
                wrote_csv = TRUE;
            } else if (s->output_dir && matches(buf, "output directory")) {
                add_to_text(&new_sav_txt, &len,
                    "output directory = %s\r\n", s->output_dir);
                wrote_output = TRUE;
            } else if (matches(buf, "next request number")) {
                add_to_text(&new_sav_txt, &len,
                    "next request number = %d\r\n", s->req_num);
                wrote_req_num = TRUE;
            } else if (matches(buf, "next request id")) {
                add_to_text(&new_sav_txt, &len,
                    "next request id = %d\r\n", s->req_id);
                wrote_req_id = TRUE;
            } else if (s->station_code && matches(buf, "station code")) {
                add_to_text(&new_sav_txt, &len,
                    "station code = %s\r\n", s->station_code);
                wrote_station_code = TRUE;
            } else {
                add_to_text(&new_sav_txt, &len, "%s", buf);
            }
        }
        FCLOSE(fp);

        if (s->csv_dir && !wrote_csv)
            add_to_text(&new_sav_txt, &len,
                "csv directory = %s\r\n", s->csv_dir);
        if (s->output_dir && !wrote_output)
            add_to_text(&new_sav_txt, &len,
                "output directory = %s\r\n", s->output_dir);
        if (!wrote_req_num && s->req_num > 1)
            add_to_text(&new_sav_txt, &len,
                "next request number = %d\r\n", s->req_num);
        if (!wrote_req_id && s->req_id > 1)
            add_to_text(&new_sav_txt, &len,
                "next request id = %d\r\n", s->req_id);
        if (s->station_code && !wrote_station_code)
            add_to_text(&new_sav_txt, &len,
                "station code = %s\r\n", s->station_code);

        fp = FOPEN(sav_file, "w");
        fprintf(fp, "%s", new_sav_txt);
        FCLOSE(fp);
    } else {
        sav_file = MALLOC(sizeof(char)*(strlen(get_asf_share_dir())+32));
        sprintf(sav_file, "%s/req_settings.txt", get_asf_share_dir());
        // write a new file
        FILE *fp = FOPEN(sav_file, "w");
        if (s->csv_dir)
            fprintf(fp, "csv directory = %s\r\n", s->csv_dir);
        if (s->output_dir)
            fprintf(fp, "output directory = %s\r\n", s->output_dir);
        fprintf(fp, "next request number = %d\r\n", s->req_num);
        fprintf(fp, "next request id = %d\r\n", s->req_id);
        if (s->station_code)
            fprintf(fp, "station code = %s\r\n", s->station_code);
        FCLOSE(fp);
    }
    free(sav_file);
}

static char *mkstr(const char *s)
{
    char *ret=NULL;
    if (strlen(s) > 0)
        ret = STRDUP(s);
    return ret;
}

Settings *settings_new_from_gui()
{
    Settings *s = settings_new();

    s->csv_dir = settings_get_csv_dir();
    s->output_dir = settings_get_output_dir();
    s->req_num = atoi(get_string_from_entry("next_request_number_entry"));
    s->req_id = atoi(get_string_from_entry("next_request_id_entry"));
    s->station_code = STRDUP(settings_get_station_code());

    return s;
}

int save_settings()
{
    Settings *s = settings_new_from_gui();
    if (!fileExists(s->csv_dir)) {
        message_box("CSV Directory doesn't exist.");
        return FALSE;
    }
    if (!fileExists(s->output_dir)) {
        message_box("Output Directory doesn't exist.");
        return FALSE;
    }
    settings_save(s);
    settings_free(s);
    return TRUE;
}

char *settings_get_csv_dir()
{
    return mkstr(get_string_from_entry("csv_directory_entry"));
}

char *settings_get_output_dir()
{
    return mkstr(get_string_from_entry("output_directory_entry"));
}

static void set_settings_saved_label(const char *txt)
{
    GtkWidget *w = get_widget_checked("settings_saved_label");
    gtk_label_set_text(GTK_LABEL(w), txt);
}

SIGNAL_CALLBACK void on_change_page(GtkWidget *widget)
{
    set_settings_saved_label("");
    put_string_to_label("generate_label", "");
    update_output_file();
}

SIGNAL_CALLBACK void on_save_button_clicked(GtkWidget *widget)
{
    if (save_settings()) {
        set_settings_saved_label("Settings Saved.");
        populate_csvs();
    }
}

int settings_get_next_req_id(void)
{
    Settings *s = settings_load();
    int id = s->req_id;
    settings_free(s);
    return id;
}

void settings_set_next_req_id_and_incr_req_num(int req_id)
{
    if (settings_get_request_type() != ON_DEMAND_LEVEL_0) {
        Settings *s = settings_load();
        if (req_id <= s->req_id)
            printf("*** New request id is smaller!?\n");
        s->req_id = req_id;
        ++s->req_num;
        settings_save(s);
        apply_settings_to_gui(s);
        settings_free(s);
        update_output_file();
    }
}

int settings_get_is_emergency()
{
    GtkWidget *w = get_widget_checked("emergency_checkbutton");
    return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));
}

int settings_get_sequence_number()
{
    return 1;
}

void settings_set_request_type(int request_type)
{
    GtkLabel *l = GTK_LABEL(get_widget_checked("request_type_label"));

    if (request_type==OBSERVATION_REQUEST)
        gtk_label_set_text(l, "Observation Request");
    else if (request_type==ACQUISITION_REQUEST)
        gtk_label_set_text(l, "Acquisition Request");
    else if (request_type==ON_DEMAND_LEVEL_0)
        gtk_label_set_text(l, "On-Demand Level 0");
    else if (request_type==UNSELECTED_REQUEST_TYPE)
        gtk_label_set_text(l, "???");
    else
        assert(0);
}

int settings_get_request_type()
{
    GtkLabel *l = GTK_LABEL(get_widget_checked("request_type_label"));
    const char *s = gtk_label_get_text(l);
    if (strcmp(s, "Observation Request") == 0)
        return OBSERVATION_REQUEST;
    else if (strcmp(s, "Acquisition Request") == 0)
        return ACQUISITION_REQUEST;
    else if (strcmp(s, "On-Demand Level 0") == 0)
        return ON_DEMAND_LEVEL_0;
    return UNSELECTED_REQUEST_TYPE;
}

const char *settings_get_station_code()
{
    static char *station_code=NULL;
    if (!station_code) {
        Settings *s = settings_load();
        station_code = STRDUP(s->station_code);
        settings_free(s);
    }
    return station_code;
}

long settings_get_start_date()
{
    char *s = get_string_from_entry("start_date_entry");
    if (strlen(s) > 0) {
        return atol(s);
    } else {
        return -1;
    }
}

void settings_set_start_date(long l)
{
    block_processing=TRUE;
    if (l < 0)
        put_string_to_entry("start_date_entry", "");
    else
        put_long_to_entry("start_date_entry", l);
    block_processing=FALSE;
}

long settings_get_end_date()
{
    char *s = get_string_from_entry("end_date_entry");
    if (strlen(s) > 0) {
        return atol(s);
    } else {
        return -1;
    }
}

void settings_set_end_date(long l)
{
    block_processing=TRUE;
    if (l < 0)
        put_string_to_entry("end_date_entry", "");
    else
        put_long_to_entry("end_date_entry", l);
    block_processing=FALSE;
}

SIGNAL_CALLBACK void on_emergency_checkbutton_toggled(GtkWidget *w)
{
    update_output_file();
    gui_process(FALSE);
}

SIGNAL_CALLBACK void on_start_date_entry_changed(GtkWidget *w)
{
    update_output_file();
    gui_process(FALSE);
}

SIGNAL_CALLBACK void on_end_date_entry_changed(GtkWidget *w)
{
    update_output_file();
    gui_process(FALSE);
}
