
#include <algorithm>  // for clamp
#include <cmath>      // for cos, sin, fmod, acos, asin, sqrt, atan2, pow
#include <cstdint>    // for int32_t
#include <exception>  // for exception
#include <functional> // for placeholders
#include <regex>      // for match_results<>::_Base_type
#include <string.h>   // for memcpy
#include <string>     // for string, allocator, operator+, to_string, char_traits
#include <time.h>     // for tm, timespec, localtime
#include <vector>     // for vector
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstdlib>
using namespace std;


#define PI 3.14159265
#define light 299792458.
#define one_over_c 0.0033356
#define R2D 180. / PI
#define D2R PI / 180.
#define TAU 2 * PI

//Defining some CHIME related values
#define _num_elements 2048
#define _num_beams 1
#define _feed_sep_NS 0.3048 
#define _feed_sep_EW 22.0
#define _inst_lat 49.32070922
#define _inst_long -119.62367743

struct beamCoord { 
    float ra[10];
    float dec[10];
    float scaling[10]; // specified the size
};


void calculate_phase(const beamCoord& beam_coord, timespec time_now, float freq_now, float* gains, float* output) 
{
    printf("Function started!\n");
    ofstream myfile2;
    myfile2.open ("feed_pos.txt");
    float FREQ = freq_now;
    struct tm* timeinfo;
    timeinfo = gmtime(&time_now.tv_sec);    //Get UTC from unix time
    printf ("Time when beamforming: %s", asctime(timeinfo)); 

    //Calculating JD from the UTC time
    uint year = timeinfo->tm_year + 1900;
    uint month = timeinfo->tm_mon + 1;
    if (month < 3) {
        month = month + 12;
        year = year - 1;
    }
    uint day = timeinfo->tm_mday;
    uint hour = timeinfo->tm_hour;
    //printf("Year=%d, month=%d, day=%d Hour=%d Min=%d Sec=%d\n", year, month, day, hour, timeinfo->tm_min, timeinfo->tm_sec);
    float JD = 2 - int(year / 100.) + int(int(year / 100.) / 4.) + int(365.25 * year)
               + int(30.6001 * (month + 1)) + day + 1720994.5;

    //Calculating GMST and thus LST (LST = GMST + longitude)
    double T = (JD - 2451545.0)
               / 36525.0; // Works if time after year 2000, otherwise T is -ve and might break
    double T0 = fmod((6.697374558 + (2400.051336 * T) + (0.000025862 * T * T)), 24.); //This takes into account the preceesion of equinox since J2000 
    double UT = (hour) + (timeinfo->tm_min / 60.)
               + (timeinfo->tm_sec + time_now.tv_nsec / 1.e9) / 3600.;
    double GST = fmod((T0 + UT * 1.002737909), 24.);
    double LST = GST + _inst_long / 15.;
    while (LST < 0) {
        LST = LST + 24;
    }
    LST = fmod(LST, 24); 
    printf("JD = %f\n", JD);
    printf("T = %lf\n", T);
    printf("T0 = %lf\n", T0);
    printf("UT = %lf\n", UT);
    printf("GST = %lf\n", GST);
    printf("LST = %lf\n", LST*15.);
    for (int b = 0; b < _num_beams; b++) {
        // If scaling is 1, then phases = gains. If scaling not equal to 1, then whatever the scaling value, doesn't 
        // matter because it won't be taken into account.
        if (beam_coord.scaling[b] == 1) {
            for (uint32_t i = 0; i < _num_elements * 2; ++i) {
                output[b * _num_elements * 2 + i] = gains[b * _num_elements * 2 + i];
            }
            continue;
        }
        double hour_angle = LST * 15. - beam_coord.ra[b];
        printf("HA = %f\n", hour_angle); 

        //Calculating alt and az from the HA
        double alt = sin(beam_coord.dec[b] * D2R) * sin(_inst_lat * D2R)
                     + cos(beam_coord.dec[b] * D2R) * cos(_inst_lat * D2R) * cos(hour_angle * D2R);
        alt = asin(std::clamp(alt, -1.0, 1.0));
        double az = (sin(beam_coord.dec[b] * D2R) - sin(alt) * sin(_inst_lat * D2R))
                    / (cos(alt) * cos(_inst_lat * D2R));
        az = acos(std::clamp(az, -1.0, 1.0));
        if (sin(hour_angle * D2R) >= 0) {
            az = TAU - az;
        }
        printf("alt = %f, az = %f\n", alt, az);
        double projection_angle, effective_angle, offset_distance;
        int cyl = 0;
        for (int i = 0; i < 4; i++) {       // loop 4 cylinders
            cyl = cyl + 1;
            myfile2 << "Cylinder ";
            myfile2 << cyl;
            myfile2 << "\n";
            for (int j = 0; j < 256; j++) { // loop 256 feeds
                //Getting feed positon for each antenna
                float dist_y = (-127.5 + j) * _feed_sep_NS; 
                float dist_x = (-1.5 + i) * _feed_sep_EW;
                myfile2 << dist_x;
                myfile2 << ", ";
                myfile2 << dist_y; 
                myfile2 << "\n";
                //Getting the real and complex phases
                projection_angle = 90 * D2R - atan2(dist_y, dist_x);
                offset_distance = sqrt(pow(dist_y, 2) + pow(dist_x, 2));
                effective_angle = projection_angle - az;
                float delay_real = cos(TAU * cos(effective_angle) * cos(alt) * offset_distance
                                       * FREQ * one_over_c);
                float delay_imag = -sin(TAU * cos(effective_angle) * cos(-alt) * offset_distance
                                        * FREQ * one_over_c);
                for (int p = 0; p < 2; p++) { // loop 2 pol
                    uint elem_id = p * 1024 + i * 256 + j;
                    // Not scrembled, assume reordering kernel has been run
                    output[(b * _num_elements + elem_id) * 2] =
                        delay_real * gains[(b * _num_elements + elem_id) * 2]
                        - delay_imag * gains[(b * _num_elements + elem_id) * 2 + 1];
                    output[(b * _num_elements + elem_id) * 2 + 1] =
                        delay_real * gains[(b * _num_elements + elem_id) * 2 + 1]
                        + delay_imag * gains[(b * _num_elements + elem_id) * 2];
                }
            }
        }
    }
    myfile2.close();
    printf("Function done!\n");
}

int main()
{

    //The array which stores the gains, initialising it with 1 + 0i for each input for now.    
    int gain_len = 2 * 2048 * _num_beams * sizeof(float);
    float * host_gain = (float*)malloc(gain_len); 
    if (host_gain == nullptr)
        throw std::runtime_error("Could not allocate memory for gains");
    int index = 0;
    for (uint b = 0; b < _num_beams * _num_elements; b++) {
        host_gain[index++] = 1;
        host_gain[index++] = 0;
        }
    //The array which will store the generated phases. Each two pair of array element stores the real and img part of the phase generated. 
    //Initializing the array elements to be 0 for now. 
    int phase_frame_len = _num_elements * _num_beams * 2 * sizeof(float);
    int scaling_frame_len = _num_beams * sizeof(float);
    float * host_phase_0 = (float*)malloc(phase_frame_len + scaling_frame_len);
    if (host_phase_0 == nullptr)
        throw std::runtime_error("Could not allocate memory for phases");
    float * host_scaling_0 = host_phase_0 + _num_elements * _num_beams * 2;
    index = 0;
    for (uint b = 0; b < _num_beams * _num_elements; b++) {
        host_phase_0[index++] = 0;
        host_phase_0[index++] = 0;
        }
    //Passing the unix time = 1613619401.449999300 for now. Time needs to be input manually.
    struct timespec time_now;
    time_now.tv_sec = 1613619401;
    time_now.tv_nsec = 449999300;

    //File to write the generated phases.
    ofstream myfile;
    myfile.open ("geo_phases.txt");

    //The source coordinates passed by the TBS. These are the coordinates wrt the precessed equinox of the time of observation. Use the j2000tonow()
    //function to compute these from the J2000 coords. Needs to be input manually for now.
    beamCoord beam_coord = {{83.95117280768928, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                            {22.027201530826435, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                            {1.5, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

    //Definig the input frequency
    float freq_now = 799.609375;
    int count = 0;
    calculate_phase(beam_coord, time_now, freq_now, host_gain, host_phase_0);
    for(index = 0; index < _num_elements * _num_beams * 2; index++){
        myfile << host_phase_0[index];
        myfile << ", ";
        count = count + 1;
    }

    printf("Phases for first 10 feeds:\n");
    for(index = 0; index <20; index++){
        printf("%f ", host_phase_0[index]);
    }
    printf("\n");
    myfile.close();
    free(host_gain);
    free(host_phase_0); 
}