#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define letters 500
#define pi 3.141592654

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

#define n_mags 29
#define n_star_mags 16
#define hexadec_places 8

#define bad_val -98.0

int n_sed;
char **sed_names;
double *sed_data,*teff,*metallicity,*ebv,*logg;

int valid_dex[n_mags];
double star_color[n_mags];

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


long long int twos_complement(long long int ii){


    int *binary;
    binary=new int[4*hexadec_places];
    int binary_places=4*hexadec_places;
    long long int remainder=ii*(-1);
    long long int local_term;
    long long int denom;
    int dex;
    int val;
    int i;
    local_term=1;
    for(i=0;i<binary_places;i++){
        local_term*=2;
    }
    for(dex=binary_places-1;dex>=0;dex--){
        local_term/=2;
        val=remainder/local_term;
        if(val>1){
            printf("Failure with binary %lld %lld %d\n",remainder,local_term,val);
            exit(1);
        }
        binary[dex]=val;
        remainder-=val*local_term;
    }


    for(i=0;i<binary_places;i++){
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
       for(i=1;i<binary_places;i++){
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
    for(dex=0;dex<binary_places;dex++){
        if(dex>0){
            local_term*=2;
        }
        remainder+=binary[dex]*local_term;
    }

    delete [] binary;

    return remainder;

}


void convert_to_hexadecimal(long long int ii, int *output){

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


void get_mag_array(long long int flag, double *mag_in, double *mag_out){
    /*
    flag is the integer flag from the catalog

    mag_in is the magnitude read in for the star

    mag_out is mag_in re-mapped to the full set of magnitudes
    that exist for each sed
    */

    int i;
    for(i=0;i<n_mags;i++){
        mag_out[i]=bad_val-10.0;
    }

    //default (sdss)
    mag_out[_2mass_j_dex]=mag_in[_star_J_dex];
    mag_out[_2mass_ks_dex]=mag_in[_star_K_dex];
    mag_out[_2mass_h_dex] = mag_in[_star_H_dex];
    mag_out[_johnson_B_dex] = mag_in[_star_B_dex];
    mag_out[_johnson_V_dex] = mag_in[_star_V_dex];
    mag_out[_sdss_u_dex] = mag_in[_star_u_dex];
    mag_out[_sdss_g_dex] = mag_in[_star_g_dex];
    mag_out[_sdss_r_dex] = mag_in[_star_r_dex];
    mag_out[_sdss_i_dex] = mag_in[_star_i_dex];
    mag_out[_sdss_z_dex] = mag_in[_star_z_dex];
    mag_out[_wise_1_dex] = mag_in[_star_w1_dex];
    mag_out[_wise_2_dex] = mag_in[_star_w2_dex];
    mag_out[_wise_3_dex] = mag_in[_star_w3_dex];
    mag_out[_wise_4_dex] = mag_in[_star_w4_dex];

    int hexadec_bits[hexadec_places];

    convert_to_hexadecimal(flag, hexadec_bits);
    int psrc=hexadec_bits[4];

    if(psrc==0 || psrc==1 || psrc==11){
        //B4 or ps1
        mag_out[_sdss_g_dex] = bad_val-10.0;
        mag_out[_sdss_r_dex] = bad_val-10.0;
        mag_out[_sdss_i_dex] = bad_val-10.0;
        mag_out[_sdss_z_dex] = bad_val-10.0;
        mag_out[_ps_g_dex] = mag_in[_star_g_dex];
        mag_out[_ps_r_dex] = mag_in[_star_r_dex];
        mag_out[_ps_i_dex] = mag_in[_star_i_dex];
        mag_out[_ps_z_dex] = mag_in[_star_z_dex];
        mag_out[_ps_y_dex] = mag_in[_star_y_dex];
    }
    else if(psrc==6 || psrc==7){
        //tycho or hipparcos
        mag_out[_johnson_B_dex] = bad_val-10.0;
        mag_out[_johnson_V_dex] = bad_val-10.0;
        mag_out[_tycho_B_dex] = mag_in[_star_B_dex];
        mag_out[_tycho_V_dex] = mag_in[_star_V_dex];
    }
    else if(psrc>12 || psrc<0){
        printf("psrc %d\n",psrc);
        printf("flag %lld\n",flag);
        exit(1);
    }


}


int fit_star_mags(double *star_mags, double *best_offset, double *err_out){
    /*
    Read in a an array of mapped star mags.
    Return the row index of the best-fitting sed, e(b-v) combination
    */


    int ii,j;
    int dex_best;
    double err,err_best;

    int n_valid=0;
    for(ii=0;ii<n_mags;ii++){
        if(star_mags[ii]>bad_val){
            n_valid++;
        }
    }


    j=0;
    for(ii=0;ii<n_mags;ii++){
        if(star_mags[ii]>bad_val){
            valid_dex[j]=ii;
            j++;
        }
    }

    double n_good=double(n_valid);
    double sed_color;

    for(j=0;j<n_valid-1;j++){
        star_color[j]=star_mags[valid_dex[j]]-star_mags[valid_dex[j+1]];
    }

    for(ii=0;ii<n_sed;ii++){
        err=0.0;
        for(j=0;j<n_valid-1;j++){
            sed_color=sed_data[ii*n_mags+valid_dex[j]]-sed_data[ii*n_mags+valid_dex[j+1]];
            err+=(sed_color-star_color[j])*(sed_color-star_color[j]);
            if(ii>0 && err>err_best){
                break;
            }
        }

        if(ii==0 || err<err_best){
            dex_best=ii;
            err_best=err;
        }
    }

    best_offset[0]=0.0;
    for(j=0;j<n_valid;j++){
        best_offset[0]+=(star_mags[valid_dex[j]]-sed_data[dex_best*n_mags+valid_dex[j]]);
    }
    best_offset[0]=best_offset[0]/n_good;

    double xx;
    err_out[0]=0.0;
    for(j=0;j<n_valid;j++){
        xx=star_mags[valid_dex[j]]-sed_data[dex_best*n_mags+valid_dex[j]]-best_offset[0];
        err_out[0]+=xx*xx;
    }

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
        printf("args: ebv_mag_grid raw_sdss_mag list_of_inputs\n");
        exit(1);
    }

    char grid_name[letters];

    char raw_grid_name[letters];

    char list_of_inputs[letters];
    int i,j;

    for(i=0;i<letters-1 && argv[1][i]!=0;i++){
        grid_name[i]=argv[1][i];
    }
    grid_name[i]=0;

    for(i=0;i<letters-1 && argv[2][i]!=0;i++){
        raw_grid_name[i]=argv[2][i];
    }
    raw_grid_name[i]=0;

    for(i=0;i<letters-1 && argv[3][i]!=0;i++){
        list_of_inputs[i]=argv[3][i];
    }
    list_of_inputs[i]=0;

    printf("%s\n%s\n%s\n",grid_name,raw_grid_name,list_of_inputs);

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

    input=fopen(list_of_inputs, "r");
    for(i=0;i<n_input_files;i++){
        fscanf(input,"%s",input_files[i]);
    }
    fclose(input);


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

    input=fopen(grid_name,"r");
    for(i=0;i<n_mags+5;i++){
        fscanf(input,"%s",word);
    }
    for(i=0;fscanf(input,"%s",word)>0;i++){
        for(j=0;j<n_mags+4;j++){
            fscanf(input,"%le",&xx);
        }
    }
    fclose(input);
    n_sed=i;
    printf("n_sed %d\n",n_sed);

    sed_names=new char*[n_sed];
    sed_data=new double[n_sed*n_mags];
    teff=new double[n_sed];
    ebv=new double[n_sed];
    metallicity=new double[n_sed];
    logg=new double[n_sed];

    for(i=0;i<n_sed;i++){
        sed_names[i]=new char[letters];
    }

    input=fopen(grid_name, "r");
    for(i=0;i<n_mags+5;i++){
        fscanf(input,"%s",word);
    }
    for(i=0;i<n_sed;i++){
        fscanf(input,"%s",sed_names[i]);
        fscanf(input,"%le %le %le %le",
        &ebv[i],&teff[i],&metallicity[i],&logg[i]);

        for(j=0;j<n_mags;j++){
            fscanf(input,"%le",&sed_data[i*n_mags+j]);
        }
    }
    fclose(input);


    int n_raw_sed;
    char **raw_sed_names;
    double *raw_magnorm;
    double **raw_sdss_mags;

    input=fopen(raw_grid_name, "r");
    for(n_raw_sed=0;fscanf(input,"%s",word)>0;n_raw_sed++){
        for(i=0;i<6;i++){
            fscanf(input,"%le",&xx);
        }
    }
    fclose(input);
    printf("n_raw_sed %d\n",n_raw_sed);

    raw_sed_names=new char*[n_raw_sed];
    raw_magnorm=new double[n_raw_sed];
    raw_sdss_mags=new double*[n_raw_sed];
    for(i=0;i<n_raw_sed;i++){
        raw_sed_names[i]=new char[letters];
        raw_sdss_mags[i]=new double[5];
    }

    input=fopen(raw_grid_name,"r");
    for(i=0;i<n_raw_sed;i++){
        fscanf(input,"%s",raw_sed_names[i]);
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

    double *mapped_star_mags;
    mapped_star_mags=new double[n_mags];

    int i_chosen;
    double offset;

    char output_name[2*letters];
    int raw_dex;
    double flux_factor;

    double err;
    double t_start = double(time(NULL));
    int ct;

    // sfd control
    //FILE *control;
    //char control_name[2*letters];

    char buffer_name[2*letters];

    int *sed_to_raw_map;
    sed_to_raw_map = new int[n_sed];
    for(i=0;i<n_sed;i++){
        for(j=0;j<n_raw_sed;j++){
            if(char_same(raw_sed_names[j], sed_names[i])==1){
                    sed_to_raw_map[i]=j;
                    continue;
            }
        }
    }

    for(i_file=0;i_file<n_input_files;i_file++){
        sprintf(buffer_name,"%s_buffer.txt",input_files[i_file]);
        printf("should be unzipping %s to %s\n",input_files[i_file],buffer_name);
        sprintf(cmd,"gunzip -c %s | sed 's/,/ /g' | sed 's/  / /g' > %s",
        input_files[i_file],buffer_name);
        system(cmd);

        sprintf(output_name,"%s_ebv_grid_fit.txt", input_files[i_file]);

        // sfd control
        //sprintf(control_name,"%s_control.txt", input_files[i_file]);
        //control=fopen(control_name,"w");

        output=fopen(output_name, "w");
        fprintf(output,"#star_id ra dec mura mudec lon lat ");
        fprintf(output,"sed magnorm flux_factor E(B-V) Teff [Fe/H] log(g) ");
        fprintf(output,"lsst_u lsst_g lsst_r lsst_i lsst_z lsst_y ");
        fprintf(output,"sdss_u(ext) sdss_g(ext) sdss_r(ext) sdss_i(ext) sdss_z(ext) ");
        fprintf(output,"sdss_u(raw) sdss_g(raw) sdss_r(raw) sdss_i(raw) sdss_z(raw) ");
        fprintf(output,"magnitude_residual\n");
        input=fopen(buffer_name,"r");
        ct=-1;
        printf("looping at %e\n",double(time(NULL))-t_start);
        while(fscanf(input,"%lld",&star_id)>0){
            ct++;
            fscanf(input,"%le %le %le %le\n",
            &ra, &dec, &mura, &mudec);

            for(i=0;i<n_star_mags;i++){
                fscanf(input,"%le",&star_mags[i]);
            }
            fscanf(input,"%lld",&flag);

            get_mag_array(flag, star_mags, mapped_star_mags);
            i_chosen=fit_star_mags(mapped_star_mags, &offset, &err);

            raw_dex=sed_to_raw_map[i_chosen];

            if(raw_dex<0){
                printf("WARNING could not find raw match for\n");
                printf("%s\n",sed_names[i_chosen]);
                exit(1);
            }

            flux_factor=exp(-0.4*log(10.0)*offset);

            lonLatFromRaDec(ra, dec, &lon, &lat);

            fprintf(output,"%lld %.12le %.12le %.12le %.12le %.12le %.12le ",
            star_id, ra, dec, mura*1000.0, mudec*1000.0, lon, lat);

            fprintf(output,"%s %le %le %le ",
            sed_names[i_chosen],raw_magnorm[raw_dex]+offset,flux_factor,ebv[i_chosen]);

            fprintf(output,"%le %le %le ",teff[i_chosen], metallicity[i_chosen], logg[i_chosen]);

            fprintf(output,"%le %le %le %le %le %le ",
            sed_data[i_chosen*n_mags+_lsst_u_dex]+offset,sed_data[i_chosen*n_mags+_lsst_g_dex]+offset,
            sed_data[i_chosen*n_mags+_lsst_r_dex]+offset,sed_data[i_chosen*n_mags+_lsst_i_dex]+offset,
            sed_data[i_chosen*n_mags+_lsst_z_dex]+offset,sed_data[i_chosen*n_mags+_lsst_y_dex]+offset);

            fprintf(output,"%le %le %le %le %le ",
            sed_data[i_chosen*n_mags+_sdss_u_dex]+offset,sed_data[i_chosen*n_mags+_sdss_g_dex]+offset,
            sed_data[i_chosen*n_mags+_sdss_r_dex]+offset,sed_data[i_chosen*n_mags+_sdss_i_dex]+offset,
            sed_data[i_chosen*n_mags+_sdss_z_dex]+offset);

            fprintf(output,"%le %le %le %le %le ",
            raw_sdss_mags[raw_dex][0]+offset, raw_sdss_mags[raw_dex][1]+offset,
            raw_sdss_mags[raw_dex][2]+offset, raw_sdss_mags[raw_dex][3]+offset,
            raw_sdss_mags[raw_dex][4]+offset);

            fprintf(output,"%le ",err);
            fprintf(output,"\n");

            // sfd control
            //convert_to_hexadecimal(flag, hexadec_bits);
            //fprintf(control,"%lld %le ",star_id, err);
            //if(hexadec_bits[4]==0 || hexadec_bits[4]==1 || hexadec_bits[4]==11){
            //    fprintf(control,"PanStarrs ");
            //}
            //else{
            //    fprintf(control,"SDSS ");
            //}
            //fprintf(control,"%le %le %le %le %le\n",
            //star_mags[_star_u_dex], star_mags[_star_g_dex],
            //star_mags[_star_r_dex], star_mags[_star_i_dex],
            //star_mags[_star_z_dex]);

        }
        fclose(input);
        fclose(output);

        // sfd control
        //fclose(control);

        sprintf(cmd,"gzip -f %s",output_name);
        system(cmd);

        sprintf(cmd,"rm %s",buffer_name);
        system(cmd);
        printf("%s in %e\n",input_files[i_file],double(time(NULL))-t_start);

    }

    printf("that took %e\n",double(time(NULL))-t_start);

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