/*
  This file is part of LS² - the Localization Simulation Engine of FU Berlin.

  Copyright 2011-2013   Heiko Will, Marcel Kyas, Thomas Hillebrandt,
  Stefan Adler, Malte Rohde, Jonathan Gunthermann
 
  LS² is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  LS² is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with LS².  If not, see <http://www.gnu.org/licenses/>.

 */

//needed for log
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "ray_noise_em.h"

#define DISTRIBUTION 2
#define SDEV 15
#define MEAN 50
#define THRESHOLD -80
#define PI 3.14159265358979323846F
#define DEGREE 18000.0F

//#define FREQUENCY 868000000
//Normal wall for Frequency 868 MHz
#define WALLREDUCTION1 40
//Beton wall for Frequency 868 MHz
#define WALLREDUCTION2 100
//Glas wall for Frequency 868 MHz
#define WALLREDUCTION3 10
//If Frequency 2.4 GHz is used, WALLREDUCTION should be as follows
#define FREQUENCY 2400000000
//Normal wall for Frequency 2.4 GHz
//#define WALLREDUCTION1 60  
//Beton wall for Frequency 2.4 GHz
//#define WALLREDUCTION2 150
//Glas wall for Frequency 2.4 GHz
//#define WALLREDUCTION3 15
#define REFLECTIONCOEFFICIENT1 0.25F
#define REFLECTIONCOEFFICIENT2 0.25F
#define REFLECTIONCOEFFICIENT3 0.05F
#define LENGTHREDUCTION 1
#define WALLTHICKNESS 0.05F
#define WALLCROSSBLOCK 1
#define MAX_DEPTH 32
//Scale: 1000 Pixel = 50 m (Like the Computer Sience Building)
#define SCALE 3
#define TRANSMITPOWER 10



/*! The name of the wall file. */
static char const * ls2_ray_noise_walls = "wall_txt/2r_walls.txt";

static GOptionEntry ray_noise_arguments[] = {
        { "walls", 0, 0, G_OPTION_ARG_FILENAME, &ls2_ray_noise_walls,
          "files describing the walls", NULL },
        { NULL }
};


void __attribute__((__nonnull__))
ls2_add_ray_noise_option_group(GOptionContext *context)
{
     GOptionGroup *group;
     group = g_option_group_new("ray-noise",
                                "Parameters to the ray tracing error model",
                                "Parameters to the ray tracing error model",
                                NULL, NULL);
     g_option_group_add_entries(group, ray_noise_arguments);
     g_option_context_add_group(context, group);
}


#include "util/util_random.c"

VECTOR fzero;
VECTOR minus_one;
VECTOR plus_one;
VECTOR *restrict  wall_x;
VECTOR *restrict  wall_y;
float *wall_width;
int *wall_kind;
VECTOR wall_number;
float *length_array;
float *strength_array;
int *wall_array;
int printhelper;
int counter;
int counter2;
int counter3;
float rz;


//Calculate the Path Loss
static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
pathloss(VECTOR length, VECTOR *restrict ploss)
{
    float help = 300000000;
    VECTOR c = VECTOR_BROADCAST(&help);
    help = 4*PI;
    VECTOR p = VECTOR_BROADCAST(&help);
    help = FREQUENCY;
    VECTOR f = VECTOR_BROADCAST(&help);
    help = 20*log10f((p[0]*f[0]*(length[0]/SCALE))/c[0]);
    *ploss = VECTOR_BROADCAST(&help);
}

//Tags every point that is part of a wall
static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
strength_on_hit(const VECTOR length, VECTOR strength, VECTOR *restrict sleft,
                VECTOR *restrict sreflect, int wnum, VECTOR angle)
{
    float wr = 0.0F;
    float rc = 0.0F;
    float inside = wall_width[wnum]/VECTOR_SIN(angle)[0];
    VECTOR ploss;
    //For debugging
    //inside = 1;
    switch(wall_kind[wnum]){
        case 1: wr = WALLREDUCTION1*inside;
                rc = REFLECTIONCOEFFICIENT1;
                break;
        case 2: wr = WALLREDUCTION2*inside;
                rc = REFLECTIONCOEFFICIENT2;
                break;
        case 3: wr = WALLREDUCTION3*inside;
                rc = REFLECTIONCOEFFICIENT3;
                break;
        default: printf("default in strength on hit \n"); 
    }
    pathloss(length, &ploss);
    strength -= ploss;
    float str = 0.1F*strength[0];
    str = powf(10,str);
    float ref = str;
    ref *= rc;
    str -= ref;
    ref = log10f(ref)*10;
    str = log10f(str)*10-wr;
    *sleft = VECTOR_BROADCAST(&str) + ploss;
    *sreflect = VECTOR_BROADCAST(&ref) + ploss;
   
    //For debugging
    /*counter3++;
    if(counter3%50000==0){
        printf("Dämpfung durch Wand kind %i ist %f \n", wall_kind[wnum], wr);
    }*/
    
}

//Get the right adress for length und strength array
static inline int
__attribute__((__always_inline__,__gnu_inline__,__const__,__artificial__))
get(int x, int y, int a)
{
    int res = (x + SIZE*y)+a*SIZE*SIZE;
    return res;
}

//tags every Point by rewriting the length and strength arrays if neccesary
static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
tag(const VECTOR ax, const VECTOR ay, const VECTOR dx, const VECTOR dy,
    const VECTOR cx __attribute__((unused)),
    const VECTOR cy __attribute__((unused)),
    int status __attribute__((__unused__)),
    const VECTOR dist,
    VECTOR *restrict length, VECTOR *restrict strength2, const VECTOR a)
{
    //counter2: how many rays are created, print in the end
    counter2++; 
    int x;
    int y;
    int b = (int) a[0];
    float r = 0.0F;
    VECTOR strength = *strength2;
    VECTOR ploss;
    while(r < dist[0]) {
        *length += plus_one;
        r++;
        x = (int) (ax[0] + r*dx[0]);
        y = (int) (ay[0] + r*dy[0]);
        //if ray leaves Area
        if(x<0||x>999||y<0||y>999) break;
        pathloss((*length), &ploss);
        strength = (*strength2) - ploss;
        //If strength is too low, break.
        if(strength[0]<THRESHOLD) {
            *strength2 = strength;
            break;
        }
        //rewrite if strength of ray is greater than saved strength and if length is smaller
        if(strength_array[get(x,y,b)] < strength[0]) {
        //Or rewrite if length is smaller, ignore strength for rewrite
        //if(length_array[get(x,y,b)] >= (*length)[0]) {
            length_array[get(x,y,b)] = (*length)[0];
            strength_array[get(x,y,b)] = strength[0];
        }
    }
}

//calculate the direction of the reflected ray
static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
ray_reflect(const VECTOR dx, const VECTOR dy, VECTOR *restrict dcx,
            VECTOR *restrict dcy, VECTOR *restrict angle)
{
    //For calculating the angle
    VECTOR alpha = dx*(*dcx)+dy*(*dcy);
    VECTOR s1 = VECTOR_SQRT(dx*dx+dy*dy);
    VECTOR s2 = VECTOR_SQRT((*dcx)*(*dcx)+(*dcy)*(*dcy));
    s1 *= s2;
    alpha /= s1;
    float beta = acosf(alpha[0]);
    if(beta>(PI/2)) beta = (beta-PI)*(-1);
    *angle = VECTOR_BROADCAST(&beta);
    //For the actual reflection
    //the perpendicular is calculated 
    VECTOR nx = (*dcy); 
    VECTOR ny = -(*dcx);
    //and normalize it 
    VECTOR norm;
    norm = VECTOR_SQRT(nx*nx+ny*ny);
    nx = nx/norm;
    ny = ny/norm;
    //Skalar product
    VECTOR skalar = dx*nx+dy*ny;
    VECTOR minus_two = minus_one + minus_one;
    *dcx = minus_two*skalar*nx + dx;
    *dcy = minus_two*skalar*ny + dy;
}

//Checks if and where a ray hits a wall
static inline int
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
ray_wall(const VECTOR ax, const VECTOR ay, const VECTOR dx, const VECTOR dy,
         VECTOR *restrict cx1, VECTOR *restrict cy1, VECTOR *restrict dcx1,
         VECTOR *restrict dcy1, VECTOR *restrict dist)
{
    VECTOR r;
    VECTOR s;
    VECTOR cx = *cx1;
    VECTOR cy = *cy1;
    VECTOR dcx = *dcx1;
    VECTOR dcy = *dcy1;
    //we need to solve the simultaneous equations
    //Ray is vertical
    if (dx[0]==0) {
        //Both are vertical
        if(dcx[0]==0){
            //Not in the same vertical line 
            if (ax[0]!=cx[0]) return 0; 
            //Find first Wallhit
            r = (cy-ay)/dy;
            if(dcy[0]<0) r = (cy-ay+dcy)/dy;
            if (r[0]<0) return 0;
            //Wallhit, Ray terminates at cx, cy
            *cy1 = ay+r*dy;
            *dist = r;
            return 2;                    
        }
        s = (ax-cx)/dcx;
        //if Wall is not hit
        if (s[0]<0||s[0]>1) return 0;
        r = (cy - ay + dcy*s)/dy;
        if (r[0]<0) return 0;
        *cx1 = ax;
        *cy1 = ay + r * dy;
        *dist = r;
        return 1;
    }
    //If ray is horizontal
    if (dy[0]==0){
        //Both are horizontal
        if(dcy[0]==0){
            //Not in the same horizontal line 
            if (ay[0]!=cy[0]) return 0; 
            //Find first Wallhit
            r = (cx-ax)/dx;
            if(dcx[0]<0) r = (cx-ax+dcx)/dx;
            if (r[0]<0) return 0;
            //Wallhit, Ray terminates at cx, cy
            *cx1 = ax+r*dx;
            *dist = r;
            return 2;                    
        }
        s = (ay-cy)/dcy;
        //Wall is not hit
        if (s[0]<0||s[0]>1) return 0;
        //where is defined by s
        r = (cx - ax + dcx*s)/dx;
        if (r[0]<0) return 0;
        *cx1 = ax + r * dx;
        *cy1 = ay;
        *dist = r;
        return 1;
    }
    //If Wall is vertical
    if(dcx[0]==0){
        r = (cx-ax)/dx;
        if (r[0]<0) return 0;
        s = (ay-cy+dy*r)/dcy;
        if (s[0]<0||s[0]>1) return 0;
        *cx1 = ax + r * dx;
        *cy1 = ay + r * dy;
        *dist = r;
        return 1;
    }
    //If Wall is horizontal
    if(dcy[0]==0){
        r = (cy-ay)/dy;
        if (r[0]<0) return 0;
        s = (ax-cx+dx*r)/dcx;
        if (s[0]<0||s[0]>1) return 0;
        *cx1 = ax + r * dx;
        *cy1 = ay + r * dy;
        *dist = r;
        return 1;
    }
    //If Wall and Ray are parallel
    if((dx[0]*dcy[0]-dy[0]*dcx[0])==0) {
        return 0;
    }
    //the "normal" case
    s = (dx*ay-dx*cy+dy*cx-dy*ax)/(dx*dcy-dy*dcx);
    if (s[0]<0||s[0]>1) return 0;
    r = (cx - ax + dcx*s)/dx;
    if (r[0]<0) return 0;
    *cx1 = ax + r * dx;
    *cy1 = ay + r * dy;
    *dist = r;
    return 1;
}


//calls ray_walls for all walls and ray_reflection
static inline int
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
ray_cross(const VECTOR ax, const VECTOR ay, const VECTOR dx, const VECTOR dy,
          VECTOR *restrict cx2, VECTOR *restrict cy2, VECTOR *restrict dcx2,
          VECTOR *restrict dcy2, VECTOR *restrict distance2,
          int *restrict wnum, VECTOR *restrict angle)
{
    VECTOR cx;
    VECTOR cy;
    VECTOR dcx;
    VECTOR dcy;
    VECTOR dist;
    float helper = SIZE*SIZE;
    *distance2 = VECTOR_BROADCAST(&helper);
    int status = 0;
    int status2 = 0;    
    for(int i = 0; i < wall_number[0]; i += 2) {
        cx = wall_x[i];
        cy = wall_y[i];
        dcx =  wall_x[i+1] - wall_x[i];
        dcy =  wall_y[i+1] - wall_y[i];
        //round for convinience
        if((dcx[0]+rz)>0&&(dcx[0]-rz)<0) {
            float help = 0.0F;
            dcx = VECTOR_BROADCAST(&help);
        }
        if((dcy[0]+rz)>0&&(dcy[0]-rz)<0) {
            float help = 0.0F;
            dcy = VECTOR_BROADCAST(&help);
        }
        status2 = ray_wall(ax,ay,dx,dy,&cx,&cy,&dcx,&dcy,&dist);
        //if wall is hit
        if (status2 != 0) {
            //and wall is closer to point
            if(dist[0]<(*distance2)[0]) {
                *distance2 = dist;
                status = status2;
                *cx2 = cx;
                *cy2 = cy;
                ray_reflect(dx,dy,&dcx,&dcy, &(*angle));
                *dcx2 = dcx;
                *dcy2 = dcy;
                *wnum = i;
            }
        } 
    }
    return status;
}

//If a Ray goes through the end of a wall, return kind of wall. This is called from ray and kills reflections.
static inline int
__attribute__((__always_inline__,__gnu_inline__,__pure__,__artificial__))
wall_cross_check(const VECTOR cx, const VECTOR cy)
{
    for(int i = 0; i < wall_number[0]; i++) {
        if (   wall_x[i][0] <= cx[0] + WALLCROSSBLOCK 
            && wall_x[i][0] >= cx[0] - WALLCROSSBLOCK
            && wall_y[i][0] <= cy[0] + WALLCROSSBLOCK
            && wall_y[i][0] >= cy[0] - WALLCROSSBLOCK) {
            return 1;
        }
    }
    return 0;
}

//The ray is handled, checking its strength, calculate the accruing rays
//Input: ARRAY of Rays, length of Array, pointer at used Ray, Pointe at next free space in Array 
static int
__attribute__((__nonnull__))
ray(VECTOR coord[],
    int max __attribute__((__unused__)),
    int depth,
    int *restrict space2)
{
    //kill Ray if not Strong enough
    if((coord[depth+5][0])<THRESHOLD) return 0;
    VECTOR cx = VECTOR_ZERO();
    VECTOR cy = VECTOR_ZERO();
    VECTOR dcx = VECTOR_ZERO();
    VECTOR dcy = VECTOR_ZERO();
    VECTOR dist = fzero;
    VECTOR angle = VECTOR_ZERO();
    //get ray from array
    VECTOR ax = coord[depth];
    VECTOR ay = coord[depth+1];
    VECTOR dx = coord[depth+2];
    VECTOR dy = coord[depth+3];
    VECTOR length = coord[depth+4];
    VECTOR strength = coord[depth+5];
    VECTOR a = coord[depth+6];
    int status = 0;
    int space = (*space2);
    int wnum = 0;
    //For debugging
    if(dx[0]==0&&dy[0]==0) {
        printf("dx = 0 und dy = 0 Bug, error beim durchaufen con coord \n");
    }
    //round for convinience
    if((dx[0] + rz) > 0 && (dx[0] - rz) <0) {
        dx = VECTOR_ZERO();
        dy = VECTOR_ONES();
    }
    if((dy[0]+rz)>0&&(dy[0]-rz)<0) {
        dx = VECTOR_ONES();
        dy = VECTOR_ZERO();
    }
    //check intersection with walls and return information. 
    //0 = no intersection, 1 = intersection, 2 = intersection with strong wall
    status = ray_cross(ax,ay,dx,dy,&cx,&cy,&dcx,&dcy,&dist,&wnum,&angle);
    //For debugging
    if(status==1&&(cx[0]>999||cy[0]>999||cx[0]<0||cy[0]<0)) {
       printf("dx = 0 Bug, wird in tag gefangen \n");
    }
    //start tag
    tag(ax,ay,dx,dy,cx,cy,status,dist,&length,&strength, a);
    //If no wall hit or strong wall hit no further rays are created 
    if(status!=1) return 0;
    //if nor==1, meaning no reflection because of egde
    int nor = 0;
    if(wall_cross_check(cx,cy)) {
        nor = 1;
    }
    //For debugging
    //nor = 1;
    VECTOR sreflect = fzero;
    VECTOR sleft = fzero;
    if(strength[0]<THRESHOLD) return 0;
    strength_on_hit(length, strength, &sleft, &sreflect, wnum, angle);
    if(sreflect[0]<THRESHOLD) nor = 1;
    VECTOR wall;
    float w = WALLTHICKNESS;
    wall = VECTOR_BROADCAST(&w);
    //kill ray if array is full
    if(space+14 > MAX_DEPTH*7) return 0;
    int res = 0;
    if(sleft[0]>=THRESHOLD) {    
        //Ray that goes through the wall
        coord[space] = cx+dx*wall;
        space++;
        coord[space] = cy+dy*wall;
        space++;
        coord[space] = dx;
        space++;
        coord[space] = dy;
        space++;
        coord[space] = length+wall;
        space++;
        coord[space] = sleft;
        space++;
        coord[space] = a;
        space++;
        res++;
    }
    if(nor == 1) {
        *space2 = space;
        return res;    
    }
    //Ray that is reflected   
    coord[space] = cx+dcx*wall;
    space++;
    coord[space] = cy+dcy*wall;
    space++;
    coord[space] = dcx;
    space++;
    coord[space] = dcy;
    space++;
    coord[space] = length+wall;
    space++;
    coord[space] = sreflect;
    space++;
    coord[space] = a;
    space++;
    res++;
    *space2 = space;
    return res;
}


//The Recursion of ray is dissolved into this function
static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__))
ray_create(const VECTOR ax, const VECTOR ay, VECTOR dx, VECTOR dy,
           VECTOR length, VECTOR strength, const VECTOR a)
{
    int go = 1;
    int depth = 0;
    int space = 7;
    //Stores new rays
    VECTOR *restrict coord;
    coord = (VECTOR*) malloc(MAX_DEPTH*7 * sizeof(VECTOR));
    if (coord == NULL) {
        //error!
        printf("ERROR malloc coord ray_create \n");
    }
    //Write ray in array
    coord[depth] = ax;
    coord[depth+1] = ay;
    coord[depth+2] = dx;
    coord[depth+3] = dy;
    coord[depth+4] = length;
    coord[depth+5] = strength;
    coord[depth+6] = a;
    
    //Realizes the Recursion
    while(go != 0) {
        if(depth > MAX_DEPTH*7) break;
        go += ray(&coord[0], MAX_DEPTH*7, depth, &space);
        depth += 7;
        go--;
    }
    free(coord);
}

//called by setup, starts n = DEGREE Rays for every anchor
static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__))
ray_start(const VECTOR ax, const VECTOR ay, VECTOR a)
{
    float helper = 0.0F;
    VECTOR length = VECTOR_BROADCAST(&helper);
    helper = TRANSMITPOWER;
    VECTOR strength = VECTOR_BROADCAST(&helper);
    VECTOR dx;
    VECTOR dy;
    for (int i = 0; i<DEGREE;i++) {
        helper = cosf((float) i * PI / (DEGREE/2.0F));
        dx = VECTOR_BROADCAST(&helper);
        helper = sinf((float) i* PI / (DEGREE/2.0F));
        dy = VECTOR_BROADCAST(&helper);
        ray_create(ax,ay,dx,dy,length,strength,a);    
    }
}


static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__))
tag_wall(void)
{
    VECTOR ax;
    VECTOR ay;
    VECTOR dx;
    VECTOR dy;
    VECTOR norm;
    int x;
    int y;
    int r = 0;
    for(int i = 0; i < wall_number[0]; i += 2) {
        r = 0;
        ax = wall_x[i];
        ay = wall_y[i];
        dx =  wall_x[i+1] - wall_x[i];
        dy =  wall_y[i+1] - wall_y[i];
        norm = VECTOR_SQRT(dx*dx+dy*dy);
        dx = dx/norm;
        dy = dy/norm;
        while(r < norm[0]+1) {
            x = (int) (ax[0] + (float) r *dx[0]);
            y = (int) (ay[0] + (float) r *dy[0]);
            wall_array[get(x,y,0)] = 10;
            r++;
        }
    }
}

//called from vector_shooter, initialize lenght and strength arrays, starts ray calculation for every anchor
static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
ray_noise_anchors(const vector2 vv[], size_t num)
{
    VECTOR anchornumber;
    VECTOR anchorx; 
    VECTOR anchory; 
    float x;
    length_array = (float *) malloc(SIZE*SIZE*num*sizeof(float *));
    if (length_array == NULL) printf("malloc length fehlgeschlagen");
    strength_array = (float *) malloc(SIZE*SIZE*num*sizeof(float *));
    if (strength_array == NULL) printf("malloc strength fehlgeschlagen");
    wall_array = (int *) malloc(SIZE*SIZE*sizeof(int *));
    if (wall_array == NULL) printf("malloc wall fehlgeschlagen");
    //initialize array
    for(size_t i = 0; i<SIZE*SIZE*num;i++) {
        length_array[i] = SIZE*SIZE;
        strength_array[i] = THRESHOLD-2;
    } 
    for(size_t i = 0; i<SIZE*SIZE;i++) {
        wall_array[i] = 0;
    }
    //Mark all Points which are part of a wall
    tag_wall();
    
    for(size_t i = 0; i<num;i++) {
        x = (float) i;
        anchornumber = VECTOR_BROADCAST(&x);
        anchorx = VECTOR_BROADCAST(&vv[i].x);
        anchory = VECTOR_BROADCAST(&vv[i].y);
        printf("%zd. Anchor at (%.0f,%.0f)\n", i+1, vv[i].x, vv[i].y);
        ray_start(anchorx, anchory, anchornumber);
    }
}

void
ray_noise_setup(const vector2 *vv __attribute__((__unused__)),
                size_t num __attribute__((__unused__)))
{
    
    printhelper = 0;
    counter = 0;
    counter2 = 0;
    counter3 = 0;
    rz = 0.0001F;
    fzero = VECTOR_BROADCASTF(0.0F);
    minus_one = VECTOR_BROADCASTF(-1.0F);
    plus_one = VECTOR_BROADCASTF(1.0F);   
    float fhelper;
    
    FILE* fp;
    fp = fopen (ls2_ray_noise_walls, "r");
    if (fp == NULL) {
        perror("Erroro opening walls file.\n");
        exit(EXIT_FAILURE);
    }
    int scan_error = 0;
    VECTOR number_of_walls = fzero;
    wall_number = fzero;
    while (fscanf(fp, "%i;", &scan_error) != EOF) {
        number_of_walls += plus_one;
    }
    rewind(fp);
    //The text may be error-prone, every wall needs 4 points and 2 values
    VECTOR vhelp;
    int now = (int) number_of_walls[0]; 
    if(now % 6 != 0) {
        now = now % 6;
        fhelper = (float) now;
        vhelp = VECTOR_BROADCAST(&fhelper);
        number_of_walls -= vhelp;
        printf("Eingabe war %i zu lang\n", now);
    }
    //set global variable wall_number
    fhelper = 6;
    vhelp = VECTOR_BROADCAST(&fhelper);
    number_of_walls /= vhelp;
    fhelper = 2;
    vhelp = VECTOR_BROADCAST(&fhelper);    
	wall_number = number_of_walls*vhelp;
	//set walls
    wall_x = (VECTOR*) malloc((size_t) wall_number[0] * sizeof(VECTOR));
    if (wall_x == NULL) printf("ERROR malloc wall_x ray_noise \n");
    
    wall_y = (VECTOR*) malloc((size_t) wall_number[0] * sizeof(VECTOR));
    if (wall_y == NULL) printf("ERROR malloc wall_y ray_noise \n");
    
    wall_width = (float *) malloc((size_t) wall_number[0] *sizeof(float *));
    if (wall_width == NULL) printf("malloc  wall_width fehlgeschlagen \n");
    
    wall_kind = (int *) malloc((size_t) wall_number[0] * sizeof(int *));
    if (wall_kind == NULL) printf("malloc wall_kind fehlgeschlagen \n");
    
    int j = 0, ihelper;
    for(int i = 0; i < number_of_walls[0]*4; i += 4) {
        scan_error = fscanf(fp, "%i;", &ihelper);
        if(scan_error == 0) printf("error scan Wand.txt");
        fhelper = (float) ihelper;
        wall_x[j] = VECTOR_BROADCAST(&fhelper);
        
        scan_error = fscanf(fp, "%i;", &ihelper);
        fhelper = (float) ihelper;
        if(scan_error == 0) printf("error scan Wand.txt");
        wall_y[j] = VECTOR_BROADCAST(&fhelper);
        
        scan_error = fscanf(fp, "%i;", &ihelper);
        fhelper = (float) ihelper;
        if(scan_error == 0) printf("error scan Wand.txt");
        wall_x[j+1] = VECTOR_BROADCAST(&fhelper);
        
        scan_error = fscanf(fp, "%i;", &ihelper);
        fhelper = (float) ihelper;
        if(scan_error == 0) printf("error scan Wand.txt");
        wall_y[j+1] = VECTOR_BROADCAST(&fhelper);
        
        scan_error = fscanf(fp, "%f;", &fhelper);
        if(scan_error == 0) printf("error scan Wand.txt");
        wall_width[j] = fhelper/1000;
        
        scan_error = fscanf(fp, "%i;", &ihelper);
        if(scan_error == 0) printf("error scan Wand.txt");
        wall_kind[j] = ihelper;
        
        j += 2;        
    }
    fclose(fp);
    //Print out Walls, for debugging    
    for(int i = 0; i <wall_number[0]; i += 2)  {
        printf("%i. Wall is (%0.0f,%0.0f) (%0.0f,%0.0f) width = %.2f kind = %i \n", 
        (i/2+1), wall_x[i][0], wall_y[i][0], wall_x[i+1][0], wall_y[i+1][0], wall_width[i], wall_kind[i]);
    }
    printf("Number of walls %.0f \n", number_of_walls[0]);
    ray_noise_anchors(vv,num);
}


static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
ray_noise_isWall(int wall[][SIZE], int zahl) {
    for (int i = 0; i<zahl; i++) {
   	    for(int j = 0; j<zahl;j++) {
       	    wall[i][j] = wall_array[i+SIZE*j];
   	    }
   	}
}



static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__(1,8)))
ray_noise_error(__m128i *restrict seed,
                const size_t anchors,
                const VECTOR *restrict distances __attribute__((__unused__)),
                const VECTOR *restrict  vx __attribute__((__unused__)),
		const VECTOR *restrict  vy __attribute__((__unused__)),
                const VECTOR tagx,
		const VECTOR tagy,
                VECTOR *restrict result)
{
    for (size_t k=0; k < anchors ; k++) {
        VECTOR rn;
        rn = gaussrand(seed, MEAN, SDEV);

        //For debbunging
        //float thelp = 0;
        //VECTOR tausend = VECTOR_BROADCAST(&thelp);
        //if no Ray hits, badluck
        //result[k] = distances[k] + rn+ tausend;
        
        int x = (int) tagx[0];
        int y = (int) tagy[0];           
        float lhelp = length_array[get(x,y, (int)k)];
        VECTOR l = VECTOR_BROADCAST(&lhelp);
        result[k] = l + rn;
    } 
    //For debugging
    if (tagx[0]==999 && tagy[0]==999 && printhelper == 0){
        
        printf("Rays created n =  %i \n", counter2);
        printf("counter3 =  %i \n", counter3);
        
        int x = 800;
        int y = 490;       
        
        printf("length anker 1 x = %i y = %i %f \n",x,y, length_array[get(x,y,0)] );
        printf("length anker 2 x = %i y = %i %f \n",x,y, length_array[get(x,y,1)] );
        printf("length anker 3 x = %i y = %i %f \n",x,y, length_array[get(x,y,2)] );
        printf("strength anker 1 x = %i y = %i %f \n",x,y, strength_array[get(x,y,0)] );
        printf("strength anker 2 x = %i y = %i %f \n",x,y, strength_array[get(x,y,1)] );
        printf("strength anker 3 x = %i y = %i %f \n",x,y, strength_array[get(x,y,2)] );
        
        //printf("return 1 =  %i \n", counter);
        //printf("who often is wall_cross_check called =  %i \n", counter3);
        printhelper++;
    }
}
