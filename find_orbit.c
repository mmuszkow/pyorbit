#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define PI         M_PI
#define POW2(x)    ((x)*(x))
#define RADIANS(x) (((x)*PI)/180.0L)
#define PIx2       (2.0L*PI)

struct satellite {
    long double a, e;
    long double rot[3];
    long double period;
    long double angle0;
    long double err;
};

struct xyzt {
    long double x, y, z, t;
};

struct xyzt* data;
int data_len;

#define POP 100
struct satellite pop[POP];
int initial_count;

long double uniform(long double a, long double b) {
    long double r01 = rand()/(long double)RAND_MAX;
    return ((b-a)*r01)+a;
}

void sat_pos(struct satellite* sat, long double t, struct xyzt* xyzt) {
    long double angle, r, ex, ey, ez, x, y, z, earth_rot;
    xyzt->t = t;
    angle = sat->angle0 + (t/sat->period) * PIx2;
    r = (sat->a * (1.0L - POW2(sat->e))) / (1.0L + sat->e * cosl(angle)); // Kepler's orbit equation
    ex = r * cosl(angle);
    ey = r * sinl(angle);
    ez = 0.0L;
    x = ex * cosl(sat->rot[1]) * cosl(sat->rot[2]) +
        ey * (cosl(sat->rot[2]) * sinl(sat->rot[0]) * sinl(sat->rot[1]) - cosl(sat->rot[0]) * sinl(sat->rot[2])) +
        ez * (cosl(sat->rot[0]) * cosl(sat->rot[2]) * sinl(sat->rot[1]) + sinl(sat->rot[0]) * sinl(sat->rot[2]));
    y = ex * cosl(sat->rot[1]) * sinl(sat->rot[2]) +
        ez * (cosl(sat->rot[0]) * sinl(sat->rot[1]) * sinl(sat->rot[2]) - cosl(sat->rot[2]) * sinl(sat->rot[0])) +
        ey * (cosl(sat->rot[0]) * cosl(sat->rot[2]) + sinl(sat->rot[0]) * sinl(sat->rot[1]) * sinl(sat->rot[2]));
    z = ez * cosl(sat->rot[0]) * cosl(sat->rot[1]) +
        ey * sinl(sat->rot[0]) * cosl(sat->rot[1]) -
        ex * sinl(sat->rot[1]);
    earth_rot = -(t/86400.0L) * PIx2;
    xyzt->x = x * cosl(earth_rot) - y * sinl(earth_rot);
    xyzt->y = y * cosl(earth_rot) + x * sinl(earth_rot);
    xyzt->z = z;
}

void compute_err(struct satellite* sat) {
    struct xyzt xyzt;
    int i;
    sat->err = 0.0L;
    for(i=0; i<data_len; i++) {
        sat_pos(sat, data[i].t, &xyzt);
        sat->err += sqrtl(POW2(data[i].x-xyzt.x) + POW2(data[i].y-xyzt.y) + POW2(data[i].z-xyzt.z));
    }
    sat->err /= (long double) data_len;
}

void randomize(struct satellite* sat) {
    sat->a = uniform(7150.0L, 7300.0L);
    sat->e = uniform(0.0L, 0.03L);
    sat->rot[0] = uniform(-PI, PI);
    sat->rot[1] = uniform(-PI, PI);
    sat->rot[2] = uniform(-PI, PI);
    sat->period = uniform(100.0L * 60.0L, 105.0L * 60.0L);
    sat->angle0 = uniform(0.0L, PIx2);
    compute_err(sat);
}

#define ANGLE_DIV 32.0L

void mutate(struct satellite* sat) {
    struct satellite mutant;
    memcpy(&mutant, sat, sizeof(struct satellite));
    switch(rand()%7) {
        case 0: mutant.a += uniform(-2.5L, 2.5L); break;
        case 1: mutant.e += uniform(-0.01L, 0.01L); break;
        case 2: {
                mutant.rot[0] += uniform(-PI/ANGLE_DIV, PI/ANGLE_DIV);
                if(mutant.rot[0] > PI) mutant.rot[0] -= PIx2;
                else if(mutant.rot[0] < -PI) mutant.rot[0] += PIx2;
                break;
                }
        case 3: {
                mutant.rot[1] += uniform(-PI/ANGLE_DIV, PI/ANGLE_DIV);
                if(mutant.rot[1] > PI) mutant.rot[1] -= PIx2;
                else if(mutant.rot[1] < -PI) mutant.rot[1] += PIx2;
                break;
                }
        case 4: {
                mutant.rot[2] += uniform(-PI/ANGLE_DIV, PI/ANGLE_DIV);
                if(mutant.rot[2] > PI) mutant.rot[2] -= PIx2;
                else if(mutant.rot[2] < -PI) mutant.rot[2] += PIx2;
                break;
                }
        case 5: mutant.period += uniform(-30.0L, 30.0L); break;
        case 6: {
                mutant.angle0 += uniform(-PI/ANGLE_DIV, PI/ANGLE_DIV);
                if(mutant.angle0 > PIx2) mutant.angle0 -= PIx2;
                else if(mutant.angle0 < 0) mutant.angle0 += PIx2;
                break;
                }
    }
    compute_err(&mutant);
    if(mutant.err < sat->err) memcpy(sat, &mutant, sizeof(struct satellite));
}

#define R_EQUAT 6378.137L
#define R_POLAR 6356.752314245L

void lon_lat2xyzt(long double lon, long double lat, long double alt, long double t, struct xyzt* xyzt) {
    long double theta, phi;
    theta = RADIANS(lon+180.0L);
    phi = RADIANS(lat+90.0L);
    xyzt->t = t;
    xyzt->x = -(R_EQUAT+alt) * sinl(phi) * cosl(theta);
    xyzt->y = -(R_EQUAT+alt) * sinl(phi) * sinl(theta);
    xyzt->z = -(R_POLAR+alt) * cosl(phi);
}

int err_func_cmp(const void* sat1, const void* sat2) {
    long double diff = ((struct satellite*) sat1)->err - ((struct satellite*) sat2)->err;
    if(diff < 0.0L) return -1;
    if(diff > 0.0L) return 1;
    return 0;
}

void find(int iters) {
    int i, j, h, m, s;
    time_t t_start, t_estim;
    
    /* initialize */
    for(i=initial_count; i<POP; i++) randomize(&pop[i]);
    for(i=1; i<POP; i++) pop[i] = pop[0];
    qsort(pop, POP, sizeof(struct satellite), err_func_cmp);

    t_start = clock();
    for(j=0; j<iters;) {
        /* mutate & sort */
        for(i=0; i<POP; i++) mutate(&pop[i]);
        qsort(pop, POP, sizeof(struct satellite), err_func_cmp);
        /* remove worst 30% */
        for(i=(int)(0.7*POP); i<POP; i++) randomize(&pop[i]);

        j++;
        /* debug */
        t_estim = (iters-j)*((clock()-t_start)/(float)j)/CLOCKS_PER_SEC;
        h = t_estim / 3600;
        m = (t_estim / 60) % 60;
        s = t_estim % 60;
        fprintf(stderr, "%d/%d err=%.21Lg ETA=%.2d:%.2d:%.2d\n", j, iters, pop[0].err, h, m, s);
    }
}

int count_lines(FILE* f) {
    char line[256];
    int count = 0;
    while(fgets(line, 255, f)) count++;
    fseek(f, 0, SEEK_SET);
    return count;
}

int main(int argc, char** argv) {
    FILE* fin;
    int i, j, id, all_count, norad_id;
    int year,month,day,hour,minu,sec;
    long double lon, lat, alt;
    struct tm tm;

    if(argc < 3) {
        printf("Usage: %s <n2yo file> <norad_id> <initial>\n", argv[0]);
        return 1;
    }
    
    /* load entries */
    fin = fopen(argv[1], "r");
    if(!fin) {
        printf("Cannot open %s\n", argv[1]);
        return 2;
    }

    norad_id = atoi(argv[2]);

    srand((unsigned int) time(NULL));

    /* count the ones with our NORAD id */
    all_count = count_lines(fin);
    data_len = 0;
    for(i = 0; i < all_count; i++) {
        fscanf(fin, "%d %d-%d-%d_%d:%d:%d %Lf %Lf %Lf\n",
                    &id, &year, &month, &day, &hour, &minu, &sec,
                    &lon, &lat, &alt);
        if(id == norad_id) data_len++;
    }
    data = (struct xyzt*) malloc(sizeof(struct xyzt) * data_len);
    fseek(fin, 0, SEEK_SET);

    /* read the data */
    for(i=0, j=0; i<all_count; i++) {
        fscanf(fin, "%d %d-%d-%d_%d:%d:%d %Lf %Lf %Lf\n",
                    &id, &year, &month, &day, &hour, &minu, &sec,
                    &lon, &lat, &alt);
        if(id != norad_id) continue;

        memset(&tm, 0, sizeof(struct tm));
        tm.tm_year = year-1900;
        tm.tm_mon = month-1;
        tm.tm_mday = day;
        tm.tm_hour = hour;
        tm.tm_min = minu;
        tm.tm_sec = sec;
        lon_lat2xyzt(lon, lat, alt, mktime(&tm), &data[j++]);
        // TODO: result is the same as in Python, but is this UTC?
        //printf("data[%d] = %Lf %Lf %Lf %Lf\n", j-1, data[j-1].x, data[j-1].y, data[j-1].z, data[j-1].t);
    }
    fclose(fin);

    /* load initial values if specified */
    initial_count = 0;
    if(argc == 4 && (fin = fopen(argv[3], "r"))) {
        all_count = count_lines(fin);
        for(i=0; i<all_count; i++) {
            fscanf(fin, "%d\t%Lf\t%Lf\t%Lf\t%Lf\t%Lf\t%Lf\t%Lf\n",
                    &id, &pop[initial_count].a, &pop[initial_count].e,
                    &pop[initial_count].rot[0], &pop[initial_count].rot[1],
                    &pop[initial_count].rot[2], &pop[initial_count].period,
                    &pop[initial_count].angle0);
            if(id == norad_id)
                compute_err(&pop[initial_count++]);

        }
        fclose(fin);
        fprintf(stderr, "%d initial values loaded\n", initial_count);
    }

    fprintf(stderr, "%d entries loaded, granularity=%.21Lg, sizeof(long double)=%lu\n",
            data_len, 1.0L/(long double)RAND_MAX, sizeof(long double));
    if(data_len > 0) {
        find(5000);
        printf("%d\t%.21Lg\t%.21Lg\t%.21Lg\t%.21Lg\t%.21Lg\t%.21Lg\t%.21Lg\n",
            norad_id, pop[0].a, pop[0].e,
            pop[0].rot[0], pop[0].rot[1], pop[0].rot[2],
            pop[0].period, pop[0].angle0);
    }

    free(data);

    return 0;
}

