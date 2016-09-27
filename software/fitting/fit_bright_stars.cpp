/*
This code reads in the outputs from ebv_independent_calc.py
and sed_mag_calc.py as well as a list of gzipped csv
files (the raw inputs provided by Dave Monet) and finds
the best fit SED, E(B-V) pair for each star in the
input csv file.  The best-fit data is written
to a file

csvName_ebv_grid_fit.txt

where csvName is the name of the csv file being fit
at any given time.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define letters 500
#define pi 3.141592654

/*
These are indices of all of the magnitudes computed
by sed_mag_calc.py
*/
#define _2mass_j_dex 0
#define _2mass_h_dex 1
#define _2mass_ks_dex 2
#define _wise_1_dex 3
#define _wise_2_dex 4
#define _wise_3_dex 5
#define _wise_4_dex 6
#define _johnson_U_dex 7
#define _johnson_B_dex 8
#define _johnson_V_dex 9
#define _hipparcos_dex 10
#define _tycho_B_dex 11
#define _tycho_V_dex 12
#define _ps_g_dex 13
#define _ps_r_dex 14
#define _ps_i_dex 15
#define _ps_z_dex 16
#define _ps_y_dex 17
#define _sdss_u_dex 18
#define _sdss_g_dex 19
#define _sdss_r_dex 20
#define _sdss_i_dex 21
#define _sdss_z_dex 22
#define _lsst_u_dex 23
#define _lsst_g_dex 24
#define _lsst_r_dex 25
#define _lsst_i_dex 26
#define _lsst_z_dex 27
#define _lsst_y_dex 28
#define _lsst_u_atm_dex 29
#define _lsst_g_atm_dex 30
#define _lsst_r_atm_dex 31
#define _lsst_i_atm_dex 32
#define _lsst_z_atm_dex 33
#define _lsst_y_atm_dex 34

/*
Each star in the input catalog only has, at most, 16 magnitudes.
They come from some cobination of the 22 bandpasses specified above
(obviously, the catalog does not contain LSST magnitudes for the stars).
In order to perform the star-to-SED fit, we read in the catalog magnitudes
of the star (indices 0 through 15), and then map those onto
a 22-element array corresponding to the 22 possible magnitudes
which have been calculated for each SED (see above; magnitudes
that do not exist for a given star are given a nonsense value).
The fit will then be calculated by comparing the color of the star
in all of the existing magnitudes to the color of each SED in those
same existing magnitudes.  The indices defined above and below
exist to standardize that process.
*/
#define _star_B_dex 0
#define _star_V_dex 1

#define _star_u_dex 2
#define _star_g_dex 3
#define _star_r_dex 4
#define _star_i_dex 5
#define _star_z_dex 6
#define _star_y_dex 7

#define _star_J_dex 8
#define _star_H_dex 9
#define _star_K_dex 10

#define _star_w1_dex 11
#define _star_w2_dex 12
#define _star_w3_dex 13
#define _star_w4_dex 14

#define _star_sst_dex 15

#define _n_allowed_colors 11

int _allowed_colors_bp1[_n_allowed_colors]={0,2,3,4,5,6,8,9,11,12,13};
double _color_min[_n_allowed_colors]={-3.0,-0.5,-0.53,-0.8,-0.32,-0.23,-1.7,-3.1,-1.03,-2.4,-1.5};
double _color_max[_n_allowed_colors]={13.3,21.0,13.5,8.3,8.2,4.6,4.3,2.6,3.1,1.0,-0.18};

#define n_mags 35
#define n_star_mags 16
#define hexadec_places 8
#define binary_places 4

#define min_mag -98.0
#define max_mag 23.0

#define KURUCZ 1
#define MLT 2
#define WD 3

#define LOG_P_KURUCZ -0.51
#define LOG_P_MLT -1.6
#define LOG_P_WD -1.6

int _n_sed;
char **_sed_names;
int *_sed_type;
double *_sed_data,*_teff,*_metallicity,*_ebv_data,*_logg;

int _valid_color_dex[_n_allowed_colors];
int _valid_mag_dex[n_star_mags];
double _star_color[_n_allowed_colors];

int _n_unq_sed=0;
int _n_ebv;
double *_unq_color,*_unq_color_max,*_unq_color_min;
int *_unq_sed_type;
char **_unq_sed_name;
int *_unq_count,*_unq_map;
int _default_mag_map[15]={8,9,18,19,20,21,22,17,0,1,2,3,4,5,6};  //from star_mag to sed mag

double power(double aa, int ee){
     double out=1.0;
     int i;
     for(i=0;i<ee;i++){
         out*=aa;
     }
     return out;
}

int char_same(char *w1, char *w2){
    int i;
    for(i=0;w1[i]!=0 && w2[i]!=0;i++){
        if(w1[i]!=w2[i]){
            return 0;
        }
    }

    if(w1[i]!=w2[i]){
        return 0;
    }

    return 1;
}

int _last_unq_dex=-1;

int get_unq_sed_dex(char *name){
    int i;
    if(_last_unq_dex>0 && _last_unq_dex<_n_unq_sed){
        if(char_same(name, _unq_sed_name[_last_unq_dex])==1){
            return _last_unq_dex;
        }
    }
    for(i=0;i<_n_unq_sed;i++){
        if(char_same(name, _unq_sed_name[i])==1){
            _last_unq_dex=i;
            return i;
        }
    }
    return -1;
}

long long int twos_complement(long long int ii){
    /*
    Compute the two's complement of a long long integer.
    This is necessary for parsing the number used to encode
    provenance information in the original catalog.
    */

    int *binary;
    binary=new int[4*hexadec_places];
    int local_binary_places=4*hexadec_places;
    long long int remainder=ii*(-1);
    long long int local_term;
    long long int denom;
    int dex;
    int val;
    int i;
    local_term=1;
    for(i=0;i<local_binary_places;i++){
        local_term*=2;
    }
    for(dex=local_binary_places-1;dex>=0;dex--){
        local_term/=2;
        val=remainder/local_term;
        if(val>1){
            printf("Failure with binary %lld %lld %d\n",remainder,local_term,val);
            exit(1);
        }
        binary[dex]=val;
        remainder-=val*local_term;
    }


    for(i=0;i<local_binary_places;i++){
        if(binary[i]==1){
            binary[i]=0;
        }
        else{
            binary[i]=1;
        }
    }

    if(binary[0]==0){
        binary[0]=1;
    }
    else{
       binary[0]=0;
       for(i=1;i<local_binary_places;i++){
           if(binary[i]==1){
               binary[i]=0;
           }
           else{
               binary[i]=1;
               break;
           }
       }
    }

    remainder=0;
    local_term=1;
    for(dex=0;dex<local_binary_places;dex++){
        if(dex>0){
            local_term*=2;
        }
        remainder+=binary[dex]*local_term;
    }

    delete [] binary;

    return remainder;

}

void convert_to_binary(int ii, int *output){
    /*
    Convert an int into an array of binary bits.  Int must be < 16
    */
    if(ii>=16){
        printf("CANNOT convert %d to binary; too big\n", ii);
        exit(1);
    }

    int denom, remainder, term;
    int i;
    denom=8;
    remainder=ii;
    for(i=binary_places-1;i>=0;i--){
        term = remainder/denom;
        if(term>1){
            printf("got a binary term that is %d\n", term);
            exit(1);
        }
        output[i] = term;
        remainder = remainder - term*denom;
        denom = denom/2;
    }
}


void convert_to_hexadecimal(long long int ii, int *output){
    /*
    Convert a long long int into an array of hexadecimal bits
    (where 10=a, 11=b, etc.)
    */

    if(ii<0){
        ii=twos_complement(ii);
    }

    long long int remainder=ii;
    long long int local_term;
    long long int denom;

    int dex;
    int val;
    int i;

    local_term=1;
    for(i=0;i<hexadec_places;i++){
        local_term*=16;
    }

    for(dex=hexadec_places-1;dex>=0;dex--){
        local_term/=16;
        val=remainder/local_term;

        if(val>=16){
            printf("WARNING cannot convert %lld to hexadecimal\n",ii);
            printf("place %d was %d\n",dex,val);
            exit(1);
        }
        output[dex]=val;
        remainder-=val*local_term;
    }
}


void get_mag_map(long long int flag, int *map_out, int *flag_out){
    int hexadec_bits[hexadec_places];

    convert_to_hexadecimal(flag, hexadec_bits);

    int i,j,k;
    int binary_bits[binary_places];
    for(i=0,j=0;i<hexadec_places;i++){
        convert_to_binary(hexadec_bits[i],binary_bits);
        for(k=0;k<binary_places;k++,j++){
            flag_out[j]=binary_bits[k];
        }
    }

    int psrc=hexadec_bits[4];

    //default (sdss)
    map_out[_star_J_dex] = _2mass_j_dex;
    map_out[_star_K_dex] = _2mass_ks_dex;
    map_out[_star_H_dex] = _2mass_h_dex;
    map_out[_star_B_dex] = _johnson_B_dex;
    map_out[_star_V_dex] = _johnson_V_dex;
    map_out[_star_u_dex] = _sdss_u_dex;
    map_out[_star_g_dex] = _sdss_g_dex;
    map_out[_star_r_dex] = _sdss_r_dex;
    map_out[_star_i_dex] = _sdss_i_dex;
    map_out[_star_z_dex] = _sdss_z_dex;
    map_out[_star_y_dex] = _ps_y_dex;
    map_out[_star_w1_dex] = _wise_1_dex;
    map_out[_star_w2_dex] = _wise_2_dex;
    map_out[_star_w3_dex] = _wise_3_dex;
    map_out[_star_w4_dex] = _wise_4_dex;

    if(psrc==0 || psrc==1 || psrc==11){
        //B4 or ps1
        map_out[_star_g_dex] = _ps_g_dex;
        map_out[_star_r_dex] = _ps_r_dex;
        map_out[_star_i_dex] = _ps_i_dex;
        map_out[_star_z_dex] = _ps_z_dex;
        map_out[_star_y_dex] = _ps_y_dex;
    }
    else if(psrc==6 || psrc==7){
        //tycho or hipparcos
        map_out[_star_B_dex] = _tycho_B_dex;
        map_out[_star_V_dex] = _tycho_V_dex;
    }
    else if(psrc>12 || psrc<0){
        printf("psrc %d\n",psrc);
        printf("flag %lld\n",flag);
        exit(1);
    }

}


int fit_star_mags(double *star_mags, int *mag_map, double *ebv_grid, double ebv_max, double *best_offset, double *err_out, int *n_valid_out){
    /*
    Read in a an array of mapped star mags.
    Return the row index of the best-fitting SED, e(b-v) combination
    */


    int ii,j,k;

    int n_valid_colors=0;
    int n_valid_mags=0;
    int ibp1,ibp2,append;
    double color;
    for(ii=0;ii<_n_allowed_colors;ii++){
        ibp1=_allowed_colors_bp1[ii];
        ibp2=ibp1+1;
        if(star_mags[ibp1]>min_mag && star_mags[ibp1]<=max_mag){
            if(star_mags[ibp2]>min_mag && star_mags[ibp2]<=max_mag){
                color=star_mags[ibp1]-star_mags[ibp2];
                if(color>=_color_min[ii] && color<=_color_max[ii]){
                    _valid_color_dex[n_valid_colors]=ii;
                    n_valid_colors++;
                    append=1;
                    for(k=0;k<n_valid_mags;k++){
                        if(_valid_mag_dex[k]==ibp1){
                            append=0;
                            break;
                        }
                    }
                    if(append==1){
                        _valid_mag_dex[n_valid_mags]=ibp1;
                        n_valid_mags++;
                    }
                    append=1;
                    for(k=0;k<n_valid_mags;k++){
                        if(_valid_mag_dex[k]==ibp2){
                            append=0;
                            break;
                        }
                    }
                    if(append==1){
                        _valid_mag_dex[n_valid_mags]=ibp2;
                        n_valid_mags++;
                    }
                }
            }
        }
    }

    n_valid_out[0]=n_valid_colors;

    if(n_valid_out[0] <= 1){
        best_offset[0]=10.0;
        err_out[0]=100.0;
        return 0;
    }

    double sed_color;

    for(j=0;j<n_valid_colors;j++){
        ibp1=_allowed_colors_bp1[_valid_color_dex[j]];
        ibp2=ibp1+1;
        _star_color[j]=star_mags[ibp1]-star_mags[ibp2];
    }

    double prior_arr[4];
    double prior;
    prior_arr[KURUCZ]= -2.0*n_valid_colors*LOG_P_KURUCZ;
    prior_arr[MLT] = -2.0*n_valid_colors*LOG_P_MLT;
    prior_arr[WD] = -2.0*n_valid_colors*LOG_P_WD;

    int i_unq;
    int i_ebv;

    int dex_best=-1;
    double err,err_best;
    double err_and_prior_best=1000000.0;

    double radius,norm,xx;
    int is_possible;

    for(i_unq=0;i_unq<_n_unq_sed;i_unq++){

        radius=-1.0;
        for(j=0;j<_n_allowed_colors;j++){
            k=i_unq*_n_allowed_colors+j;
            norm=0.5*(_unq_color_max[k]-_unq_color_min[k]);
            if(norm>radius){
                radius=norm;
            }
        }
        radius=power(radius,2);

        norm=0.0;
        for(j=0;j<n_valid_colors;j++){
            xx=_star_color[j]-_unq_color[i_unq*_n_allowed_colors+_valid_color_dex[j]];
            norm+=power(xx,2);
        }

        ii=_unq_map[i_unq*_n_ebv];
        if(_sed_type[ii]==KURUCZ){
            prior=prior_arr[KURUCZ];
        }
        else if(_sed_type[ii]==MLT){
            prior=prior_arr[MLT];
        }
        else if(_sed_type[ii]==WD){
            prior=prior_arr[WD];
        }

        is_possible=1;
        if(norm/radius>1.0 && power(sqrt(norm)-sqrt(radius),2)+prior>err_and_prior_best){
            is_possible=0;
        }

        if(dex_best>=0 && is_possible==0){
            continue;
        }

        for(i_ebv=0;i_ebv<_n_ebv;i_ebv++){
            ii=_unq_map[i_unq*_n_ebv+i_ebv];

            if(ebv_grid[ii]<=ebv_max){

                err=0.0;

                for(j=0;j<n_valid_colors;j++){
                    ibp1=_allowed_colors_bp1[_valid_color_dex[j]];
                    ibp2=ibp1+1;

                    sed_color=_sed_data[ii*n_mags+mag_map[ibp1]]-_sed_data[ii*n_mags+mag_map[ibp2]];
                    err+=(sed_color-_star_color[j])*(sed_color-_star_color[j]);
                    if(dex_best>=0 && err+prior>err_and_prior_best){
                        break;
                    }
                }

                if(dex_best<0 || err+prior<err_and_prior_best){
                    dex_best=ii;
                    err_best=err;
                    err_and_prior_best=err+prior;
                 }
            }
        }
    }

    best_offset[0]=0.0;
    for(j=0;j<n_valid_mags;j++){
        best_offset[0]+=(star_mags[_valid_mag_dex[j]]-_sed_data[dex_best*n_mags+mag_map[_valid_mag_dex[j]]]);
    }
    best_offset[0]=best_offset[0]/double(n_valid_mags);

    err_out[0]=sqrt(err_best/(n_valid_colors));

    return dex_best;

}

void lonLatFromRaDec(double ra, double dec, double *lon, double *lat){
    /*
    This is just a re-implementation of the palpay eqgal method and the
    erfa methods underlying it.
    */

    double x,y,z;
    double rr,dd;
    rr=ra*pi/180.0;
    dd=dec*pi/180.0;

    double cc=cos(dd);

    x=cos(rr)*cc;
    y=sin(rr)*cc;
    z=sin(dd);

    double rmat[3][3]={
      { -0.054875539726,-0.873437108010,-0.483834985808 },
      { +0.494109453312,-0.444829589425,+0.746982251810 },
      { -0.867666135858,-0.198076386122,+0.455983795705 }
    };

    double xp,yp,zp;
    xp=rmat[0][0]*x+rmat[0][1]*y+rmat[0][2]*z;
    yp=rmat[1][0]*x+rmat[1][1]*y+rmat[1][2]*z;
    zp=rmat[2][0]*x+rmat[2][1]*y+rmat[2][2]*z;

    double dist;
    dist=xp*xp+yp*yp;
    if(dist==0.0){
        lon[0]=0.0;
    }
    else{
        lon[0]=atan2(yp,xp);
    }

    if(zp==0.0){
        lat[0]=0.0;
    }
    else{
        lat[0]=atan2(zp,sqrt(dist));
    }

    lon[0]*=180.0/pi;
    lat[0]*=180.0/pi;

    while(lon[0]<0.0){
        lon[0]+=360.0;
    }

    while(lon[0]>360.0){
        lon[0]-=360.0;
    }

    while(lat[0]<-180.0){
        lat[0]+=360.0;
    }

    while(lat[0]>180.0){
        lat[0]-=360.0;
    }

}

int main(int iargc, char *argv[]){

    if(iargc==1){
        printf("args:\n");
        printf("-m: magnitude_grid is the output of sed_mag_calc.py\n");
        printf("-e: ebv_independent_data is the output of ebv_independent_calc.py\n");
        printf("-i: list_of_inputs is a text file listing the gzipped csv files to process\n");
        printf("-o: output_dir is the directory to put output\n");
        exit(1);
    }

    char sed_grid_name[letters];
    char raw_grid_name[letters];
    char list_of_inputs[letters];
    char output_dir[letters];
    int i,j;

    for(i=1;i<iargc;i++){
        if(argv[i][0]=='-'){
            switch(argv[i][1]){
                case 'm':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        sed_grid_name[j]=argv[i][j];
                    }
                    sed_grid_name[j]=0;
                    break;
                case 'e':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        raw_grid_name[j]=argv[i][j];
                    }
                    raw_grid_name[j]=0;
                    break;
                case 'i':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        list_of_inputs[j]=argv[i][j];
                    }
                    list_of_inputs[j]=0;
                    break;
                case 'o':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        output_dir[j]=argv[i][j];
                    }
                    output_dir[j]=0;
                    break;
            }
        }
    }

    printf("%s\n%s\n%s\n",sed_grid_name,raw_grid_name,list_of_inputs);

    FILE *input;
    char word[500];
    int n_input_files=0;
    input=fopen(list_of_inputs, "r");
    for(i=0;fscanf(input,"%s",word)>0;i++);
    fclose(input);
    n_input_files=i;


    char **input_files;
    input_files=new char*[n_input_files];
    for(i=0;i<n_input_files;i++){
        input_files[i]=new char[letters];
    }

    // read in the names of the csv files to be processed
    input=fopen(list_of_inputs, "r");
    for(i=0;i<n_input_files;i++){
        fscanf(input,"%s",input_files[i]);
    }
    fclose(input);


    // do some testing to make sure I converted to binary correctly

    printf("testing binary\n");
    int binary_bits[hexadec_places];
    convert_to_binary(8, binary_bits);
    if(binary_bits[0]!=0 || binary_bits[1]!=0 || binary_bits[2]!=0 ||
       binary_bits[3]!=1){
        printf("WARNING binary test got 8 wrong\n");
        for(i=0;i<binary_places;i++){
            printf("%d\n",binary_bits[i]);
        }
        exit(1);
    }

    convert_to_binary(12, binary_bits);
    if(binary_bits[0]!=0 || binary_bits[1]!=0 || binary_bits[2]!=1 ||
       binary_bits[3]!=1){
        printf("WARNING binary test got 12 wrong\n");
        for(i=0;i<binary_places;i++){
            printf("%d\n",binary_bits[i]);
        }
        exit(1);
    }

    convert_to_binary(7, binary_bits);
    if(binary_bits[0]!=1 || binary_bits[1]!=1 || binary_bits[2]!=1 ||
       binary_bits[3]!=0){
        printf("WARNING binary test got 7 wrong\n");
        for(i=0;i<binary_places;i++){
            printf("%d\n",binary_bits[i]);
        }
        exit(1);
    }

    // do some testing to make sure I converted
    // base 10 to base 16 properly
    printf("testing hexadec\n");

    int hexadec_bits[hexadec_places];
    long long int ii;
    ii=32;
    convert_to_hexadecimal(ii, hexadec_bits);
    for(i=0;i<hexadec_places;i++){
        if(i!=1){
            if(hexadec_bits[i]!=0){
                printf("WARNING got 32 wrong %d %d\n",i,hexadec_bits[i]);
                exit(1);
            }
        }
        else{
            if(hexadec_bits[i]!=2){
                printf("WARNING got 32 wrong %d %d\n",i,hexadec_bits[i]);
                exit(1);
            }
        }
    }


    ii=431;
    int is_wrong;
    convert_to_hexadecimal(ii,hexadec_bits);
    for(i=0;i<hexadec_places;i++){
        is_wrong=0;
        if(i==0){
            if(hexadec_bits[i]!=15){
                is_wrong=1;
            }
        }
        else if(i==1){
            if(hexadec_bits[i]!=10){
                is_wrong=1;
            }
        }
        else if(i==2){
            if(hexadec_bits[i]!=1){
                is_wrong=1;
            }
        }
        else{
            if(hexadec_bits[i]!=0){
                is_wrong=1;
            }
        }

        if(is_wrong==1){
            if(hexadec_bits[i]!=0){
                printf("WARNING got 431 wrong %d %d\n",i,hexadec_bits[i]);
                exit(1);
            }
        }
    }


    ii=218108420;
    convert_to_hexadecimal(ii,hexadec_bits);
    for(i=0;i<hexadec_places;i++){
        is_wrong=0;
        if(i==0){
            if(hexadec_bits[i]!=4){
                is_wrong=1;
            }
        }
        else if(i==2){
            if(hexadec_bits[i]!=2){
                is_wrong=1;
            }
        }
        else if(i==3){
            if(hexadec_bits[i]!=1){
                is_wrong=1;
            }
        }
        else if(i==6){
            if(hexadec_bits[i]!=13){
                is_wrong=1;
            }
        }

        if(is_wrong==1){
            printf("got %lld wrong; %d %d\n",
            ii,i,hexadec_bits[i]);
            exit(1);
        }
    }

    ii=-53;
    convert_to_hexadecimal(ii,hexadec_bits);
    for(i=0;i<hexadec_places;i++){
        is_wrong=0;

        if(i==0){
            if(hexadec_bits[i]!=11){
                is_wrong=1;
            }
        }
        else if(i==1){
            if(hexadec_bits[i]!=12){
                is_wrong=1;
            }
        }
        else{
            if(hexadec_bits[i]!=15){
                is_wrong=1;
            }
        }

        if(is_wrong==1){
            printf("got %lld wrong; %d %d\n",ii,i,hexadec_bits[i]);
        }
    }


    printf("tested all hexadec\n");


    double xx;

    // determine how many SEDs and E(B-V)s are in our grid of possible fits
    input=fopen(sed_grid_name,"r");
    for(i=0;i<n_mags+2;i++){
        fscanf(input,"%s",word);
    }
    for(i=0;fscanf(input,"%s",word)>0;i++){
        for(j=0;j<n_mags+1;j++){
            fscanf(input,"%le",&xx);
        }
    }
    fclose(input);
    _n_sed=i;
    printf("n_sed %d\n",_n_sed);

    _sed_names=new char*[_n_sed];
    _sed_type=new int[_n_sed];
    _sed_data=new double[_n_sed*n_mags];
    _ebv_data=new double[_n_sed];
    _unq_sed_name=new char*[_n_sed];

    int unq_dex;

    for(i=0;i<_n_sed;i++){
        _sed_names[i]=new char[letters];
        _unq_sed_name[i]=new char[letters];
    }

    _unq_count=new int[_n_sed];

    _n_unq_sed=0;
    input=fopen(sed_grid_name,"r");
    for(i=0;i<n_mags+2;i++){
        fscanf(input,"%s",word);
    }
    for(i=0;fscanf(input,"%s",_sed_names[i])>0;i++){

        unq_dex=get_unq_sed_dex(_sed_names[i]);
        if(unq_dex<0){
            for(j=0;j<letters;j++){
                _unq_sed_name[_n_unq_sed][j]=_sed_names[i][j];
            }
            _unq_count[_n_unq_sed]=1;
            _n_unq_sed++;
        }
        else{
            if(unq_dex>=_n_unq_sed){
                printf("WARNING unq dex too big whle making names\n");
                exit(1);
            }
            _unq_count[unq_dex]+=1;
        }

        for(j=0;j<n_mags+1;j++){
            fscanf(input,"%le",&xx);
        }
    }
    fclose(input);

    for(i=0;i<_n_unq_sed;i++){
        if(_unq_count[i]!=_unq_count[0]){
            printf("WARNING: grid not regular %d\n",_unq_count[0]);
            for(j=0;j<_n_unq_sed;j++){
                if(_unq_count[j]!=_unq_count[0]){
                    printf("%s %d\n",_unq_sed_name[j],_unq_count[j]);
                }
            }
            exit(1);
        }
    }

    printf("%d unique seds\n",_n_unq_sed);

    _unq_color=new double[_n_unq_sed*_n_allowed_colors];
    _unq_sed_type=new int[_n_unq_sed];
    double ebv_min=100.0;

    int ibp1, ibp2, i_color;
    double color;
    _unq_color_max=new double[_n_unq_sed*_n_allowed_colors];
    _unq_color_min=new double[_n_unq_sed*_n_allowed_colors];

    for(i=0;i<_n_unq_sed*_n_allowed_colors;i++){
        _unq_color_max[i]=-1000000.0;
        _unq_color_min[i]=1000000.0;
    }

    _n_ebv=_unq_count[0];
    _unq_map=new int[_n_unq_sed*_n_ebv];

    for(i=0;i<_n_unq_sed;i++){
        _unq_count[i]=0;
    }

    // read in SEDs and E(B-V)s
    input=fopen(sed_grid_name, "r");
    for(i=0;i<n_mags+2;i++){
        fscanf(input,"%s",word);
    }
    for(i=0;i<_n_sed;i++){
        fscanf(input,"%s",_sed_names[i]);
        if(_sed_names[i][0]=='k'){
            _sed_type[i]=KURUCZ;
        }
        else if(_sed_names[i][0]=='l'){
            _sed_type[i]=MLT;
        }
        else if(_sed_names[i][0]=='b'){
            _sed_type[i]=WD;
        }
        else{
           printf("Do not know type for %s (%d)\n",_sed_names[i],i);
           exit(1);
        }
        fscanf(input,"%le",&_ebv_data[i]);
        if(_ebv_data[i]<ebv_min){
            ebv_min=_ebv_data[i];
        }

        for(j=0;j<n_mags;j++){
            fscanf(input,"%le",&_sed_data[i*n_mags+j]);
        }

        unq_dex=get_unq_sed_dex(_sed_names[i]);

        if(unq_dex<0 || unq_dex>=_n_unq_sed){
            printf("WARNING unq_dex %d; but %d\n",unq_dex,_n_unq_sed);
            exit(1);
        }

        if(_unq_count[unq_dex]>=_n_ebv){
            printf("WARNING unq count exceeded _n_ebv %d %d %s %s\n",
            _n_ebv,_unq_count[unq_dex],_unq_sed_name[unq_dex],_sed_names[i]);
            exit(1);
        }

        _unq_map[unq_dex*_n_ebv+_unq_count[unq_dex]]=i;
        _unq_count[unq_dex]+=1;
        _unq_sed_type[unq_dex]=_sed_type[i];

        for(i_color=0;i_color<_n_allowed_colors;i_color++){
            ibp1=_default_mag_map[_allowed_colors_bp1[i_color]];
            ibp2=_default_mag_map[_allowed_colors_bp1[i_color]+1];
            if(ibp1>=n_mags || ibp2>=n_mags){
               printf("WARNING ibp1/2 out of bounds\n");
               exit(1);
            }
            color = _sed_data[i*n_mags+ibp1]-_sed_data[i*n_mags+ibp2];
            if(color<_unq_color_min[unq_dex*_n_allowed_colors+i_color]){
                _unq_color_min[unq_dex*_n_allowed_colors+i_color]=color;
            }
            if(color>_unq_color_max[unq_dex*_n_allowed_colors+i_color]){
                _unq_color_max[unq_dex*_n_allowed_colors+i_color]=color;
            }
        }

    }
    fclose(input);

    for(i=0;i<_n_unq_sed;i++){
        for(i_color=0;i_color<_n_allowed_colors;i_color++){
            _unq_color[i*_n_allowed_colors+i_color]= \
            0.5*(_unq_color_min[i*_n_allowed_colors+i_color]+_unq_color_max[i*_n_allowed_colors+i_color]);
        }
    }

    int n_raw_sed;
    char **raw_sed_names;
    double *raw_magnorm;
    double **raw_sdss_mags;

    // read in the raw SDSS magnitudes and magNorms of our SEDs
    input=fopen(raw_grid_name, "r");
    for(n_raw_sed=0;fscanf(input,"%s",word)>0;n_raw_sed++){
        for(i=0;i<9;i++){
            fscanf(input,"%le",&xx);
        }
    }
    fclose(input);
    printf("n_raw_sed %d\n",n_raw_sed);

    raw_sed_names=new char*[n_raw_sed];
    raw_magnorm=new double[n_raw_sed];
    raw_sdss_mags=new double*[n_raw_sed];
    _teff=new double[n_raw_sed];
    _metallicity=new double[n_raw_sed];
    _logg=new double[n_raw_sed];
    for(i=0;i<n_raw_sed;i++){
        raw_sed_names[i]=new char[letters];
        raw_sdss_mags[i]=new double[5];
    }

    input=fopen(raw_grid_name,"r");
    for(i=0;i<n_raw_sed;i++){
        fscanf(input,"%s",raw_sed_names[i]);
        fscanf(input,"%le %le %le",&_teff[i],&_metallicity[i],&_logg[i]);
        fscanf(input,"%le",&raw_magnorm[i]);
        for(j=0;j<5;j++){
            fscanf(input,"%le",&raw_sdss_mags[i][j]);
        }
    }
    fclose(input);

    int i_file;
    char cmd[2*letters];

    long long int star_id;
    double ra,dec,lon,lat;
    double mura,mudec;
    double *star_mags;
    star_mags=new double[n_star_mags];
    long long int flag;


    FILE *output;

    int *mag_map;
    mag_map = new int[n_star_mags-1];

    int i_chosen;
    double offset;

    char output_name[2*letters];
    int raw_dex;
    double flux_factor;

    double err;
    int ct;

    char buffer_name[2*letters];

    // construct a mapping so we can associate our SED, E(B-V) grid to
    // the raw SDSS magnitudes and magNorms
    int *sed_to_raw_map;
    sed_to_raw_map = new int[_n_sed];
    for(i=0;i<_n_sed;i++){
        sed_to_raw_map[i]=-1;
        for(j=0;j<n_raw_sed;j++){
            if(char_same(raw_sed_names[j], _sed_names[i])==1){
                    sed_to_raw_map[i]=j;
                    break;
            }

        }
        if(sed_to_raw_map[i]<0){
            printf("WARNING could not find raw_map for %s\n",
            _sed_names[i]);
            exit(1);
        }
    }

    char ebv_max_name[2*letters];

    double *ebv_max=NULL;
    double ebv_max_scalar;
    int n_ebv_max,i_line;
    double t_start = double(time(NULL));
    int total_ct=0;
    int n_colors;
    int monet_bits[hexadec_places*binary_places];
    char highest_name[letters];
    // Loop over the csv files, performing the fits of all stars in those files
    for(i_file=0;i_file<n_input_files;i_file++){

        for(i=0,j=0;i<letters && input_files[i_file][j]!=0;j++){
            if(input_files[i_file][j]=='/'){
                i=0;
            }
            else{
                highest_name[i]=input_files[i_file][j];
                i++;
            }
        }
        highest_name[i]=0;

        sprintf(buffer_name,"%s_buffer.txt",input_files[i_file]);

        // unzip the csv file and use sed to make it ' ' delimited, rather than ',' delimited
        printf("should be unzipping %s to %s\n",input_files[i_file],buffer_name);
        sprintf(cmd,"gunzip -c %s | sed 's/,/ /g' | sed 's/  / /g' > %s",
        input_files[i_file],buffer_name);
        system(cmd);

        sprintf(ebv_max_name,"%s_ebv_max.txt",input_files[i_file]);

        sprintf(cmd,"python generate_ebv_max.py %s %s",input_files[i_file],ebv_max_name);
        system(cmd);

        input=fopen(ebv_max_name, "r");
        for(n_ebv_max=0;fscanf(input,"%le",&xx)>0;n_ebv_max++);
        fclose(input);
        if(ebv_max!=NULL){
            delete [] ebv_max;
        }
        ebv_max=new double[n_ebv_max];
        input=fopen(ebv_max_name, "r");
        for(i=0;i<n_ebv_max;i++){
            fscanf(input,"%le",&ebv_max[i]);
        }

        sprintf(output_name,"%s/%s_ebv_grid_fit.txt", output_dir, highest_name);

        output=fopen(output_name, "w");

        // read in the ' ' delimited file created with sed
        input=fopen(buffer_name,"r");
        ct=-1;
        printf("looping at %e\n",double(time(NULL))-t_start);

        // loop over all of the stars in the ' ' delimited file
        for(i_line=0;fscanf(input,"%lld",&star_id)>0;i_line++){
            if(i_line>=n_ebv_max){
                printf("WARNING i_line %d n_ebv_max %d\n",i_line,n_ebv_max);
                exit(1);
            }
            total_ct++;
            ct++;
            if(total_ct%1000==0){
                printf("did %d in %e hours\n",total_ct,
                (double(time(NULL))-t_start)/3600.0);
            }
            fscanf(input,"%le %le %le %le\n",
            &ra, &dec, &mura, &mudec);

            for(i=0;i<n_star_mags;i++){
                fscanf(input,"%le",&star_mags[i]);
            }
            fscanf(input,"%lld",&flag);

            // map the star's magnitudes onto the SED magnitudes
            get_mag_map(flag, mag_map, monet_bits);

            // make sure that ebv_max fits in the grid of ebv we have input
            if(ebv_max[i_line]>ebv_min){
                ebv_max_scalar=ebv_max[i_line];
            }
            else{
                ebv_max_scalar=ebv_min+1.0e-6;
            }

            // choose the SED, E(B-V) pair that best matches the star's colors
            i_chosen=fit_star_mags(star_mags, mag_map, _ebv_data, ebv_max_scalar, &offset, &err, &n_colors);

            raw_dex=sed_to_raw_map[i_chosen];

            if(raw_dex<0){
                printf("WARNING could not find raw match for\n");
                printf("%s\n",_sed_names[i_chosen]);
                exit(1);
            }

            flux_factor=exp(-0.4*log(10.0)*offset);

            lonLatFromRaDec(ra, dec, &lon, &lat);

            fprintf(output,"%lld %.12le %.12le %.12le %.12le %.12le %.12le ",
            star_id, ra, dec, mura*1000.0, mudec*1000.0, lon, lat);

            fprintf(output,"%s %le %le %le ",
            _sed_names[i_chosen],raw_magnorm[raw_dex]+offset,flux_factor,_ebv_data[i_chosen]);

            fprintf(output,"%le %le %le %le %le %le ",
            _sed_data[i_chosen*n_mags+_lsst_u_dex]+offset,_sed_data[i_chosen*n_mags+_lsst_g_dex]+offset,
            _sed_data[i_chosen*n_mags+_lsst_r_dex]+offset,_sed_data[i_chosen*n_mags+_lsst_i_dex]+offset,
            _sed_data[i_chosen*n_mags+_lsst_z_dex]+offset,_sed_data[i_chosen*n_mags+_lsst_y_dex]+offset);

            fprintf(output,"%le %le %le %le %le %le ",
            _sed_data[i_chosen*n_mags+_lsst_u_atm_dex]+offset,_sed_data[i_chosen*n_mags+_lsst_g_atm_dex]+offset,
            _sed_data[i_chosen*n_mags+_lsst_r_atm_dex]+offset,_sed_data[i_chosen*n_mags+_lsst_i_atm_dex]+offset,
            _sed_data[i_chosen*n_mags+_lsst_z_atm_dex]+offset,_sed_data[i_chosen*n_mags+_lsst_y_atm_dex]+offset);

            for(i=0;i<n_star_mags;i++){
                fprintf(output,"%le ",star_mags[i]);
            }
            for(i=0;i<hexadec_places*binary_places;i++){
                fprintf(output,"%d ",monet_bits[i]);
            }
            fprintf(output,"%le %s ",err,highest_name);
            fprintf(output,"%d\n",n_colors);

        }
        fclose(input);
        fclose(output);

        sprintf(cmd,"gzip -f %s",output_name);
        system(cmd);

        sprintf(cmd,"rm %s",buffer_name);
        system(cmd);
        sprintf(cmd,"rm %s",ebv_max_name);
        system(cmd);
        printf("%s in %e\n",input_files[i_file],double(time(NULL))-t_start);

    }

    printf("that took %e to do %d\n",double(time(NULL))-t_start,total_ct);

}

/*
 char **mag_names;
    mag_names=new char*[n_star_mags];
    for(i=0;i<n_star_mags;i++){
        mag_names[i]=new char[10];
    }

    sprintf(mag_names[0],"B");
    sprintf(mag_names[1],"V");
    sprintf(mag_names[2],"u");
    sprintf(mag_names[3],"g");
    sprintf(mag_names[4],"r");
    sprintf(mag_names[5],"i");
    sprintf(mag_names[6],"z");
    sprintf(mag_names[7],"y");
    sprintf(mag_names[8],"J");
    sprintf(mag_names[9],"H");
    sprintf(mag_names[10],"K");
    sprintf(mag_names[11],"w1");
    sprintf(mag_names[12],"w2");
    sprintf(mag_names[13],"w3");
    sprintf(mag_names[14],"w4");
    sprintf(mag_names[15],"sst");

    char **cat_map;
    cat_map=new char*[13];
    for(i=0;i<13;i++){
        cat_map[i]=new char[100];
    }
    sprintf(cat_map[0],"B4-0");
    sprintf(cat_map[1],"B4-1");
    sprintf(cat_map[2],"2mass");
    sprintf(cat_map[3],"sdss");
    sprintf(cat_map[4],"ppmxl");
    sprintf(cat_map[5],"ucac4");
    sprintf(cat_map[6],"tycho");
    sprintf(cat_map[7],"hipparcos");
    sprintf(cat_map[8],"wise");
    sprintf(cat_map[9],"lepine_shara");
    sprintf(cat_map[10],"apass");
    sprintf(cat_map[11],"ps1");
    sprintf(cat_map[12],"2mass_ext");
*/
