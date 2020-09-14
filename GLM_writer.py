import json

# parameters = data['initial_parameters']


# start = parameters['STARTTIME']
# stop = parameters['STOPTIME']
#
# f.write('#define STARTTIME="{}-{}-{}T{}:{}:{}-00:00"\n'.format(start['year'],start['month'],start['day'],start['hour'],start['min'],start['sec']))
# f.write('#define STOPTIME="{}-{}-{}T{}:{}:{}-00:00"\n'.format(stop['year'],stop['month'],stop['day'],stop['hour'],stop['min'],stop['sec']))
# f.write('#define TIMEZONE="{}"\n'.format(parameters['TIMEZONE']))
# f.write('#define RECORDINTERVAL="{}"\n'.format(parameters['RECORDINTERVAL']))
# f.write('#define TMY3_FILE="{}"\n'.format(parameters['TMY3_FILE']))



def write_glm_flex(name_glm_file):
    # This function can be used on any GLM file, whereas write_recovery_ieee37 function is hardcoded for the IEEE37 GLM
    # It replaces a faulty string if it exists
    # It appends "'#include 'ESCs.glm'" to the files
    file = name_glm_file

    #read input file
    fin = open(file,'rt')

    #read file contents to string
    data = fin.read()

    # If faulty string in GLM file, replace it.
    if 'global enumeration climate_impact_zone NONE;' in data:

        #replace all occurrences of the required string
        data = data.replace('global enumeration climate_impact_zone NONE;', 'global enumeration climate_impact_zone 0;')

        #close the input file
        fin.close()

        #open the input file in write mode
        fin = open(file, "wt")

        #overwrite the input file with the resulting data
        fin.write(data)

    #close the file
    fin.close()

    # Open in append mode and append "#include 'ESCs.glm'"
    with open(file, "a") as myfile:
        myfile.write("\n")
        myfile.write('#include "ESCs.glm"')
        myfile.write("\n")

    myfile.close()

def write_recovery_ieee37(name_glm_file):
    # Write an updated recovery_ieee37.glm with updated start/stop times. Not required in GRIP plaform
    with open('config.json','r') as f:
        data = json.load(f)

    f = open(name_glm_file, "w+", newline='')

    f.write('#define STARTTIME="{}"\n'.format(data['starttime']))
    f.write('#define STOPTIME="{}"\n'.format(data['stoptime']))
    f.write('#define RECORDINTERVAL="1"\n')
    f.write('#define TMY3_FILE="{}"\n'.format(data['weather_filename']))

    f.write('''
#define climate_object="yes"

#define LOADS="on" // enabling explicit loads
#define SOLAR="on"
#set dateformat=ISO8601
//#define PVPROB //enabling a distribution of houses with and without solar
//if PVPROB is not defined assumes if SOLAR=on all houses have PV panels
//#define PVAREA=400 sf //area of each PV panel

// SETTINGS
//#set relax_naming_rules=1;
clock {
	//timezone ${TIMEZONE};
	starttime "${STARTTIME}";
	stoptime "${STOPTIME}";
};

// Import powerflow module
module powerflow {
	solver_method FBS;
};

// Import generator module for solar and inverter
module generators;

// Import climate module for solar weather data
module climate;

// IMPORT PYTHON MODULE
module ES_controller;
//module ES_controller_with_noise;

// module for recorder objects:
module tape;

// INCLUDED FILES
#include "library.glm"
#include "feeder.glm"
#include "loads.glm"
//#include "explicit_loads.glm"
#include "ESCs.glm"

#ifdef climate_object
object climate {
    name climate_obj;
    tmyfile ${TMY3_FILE};
}
#endif // climate object

object regulator {
	name reg_799781;
	to node_781;
	phases "ABC";
	from node_799;
	configuration regulator_configuration_79978101;
};

// Recorder to save outputs in CSV file (here only used to create a fixed timestep, recordings are from python module)
object recorder {
    parent line_701_702;
	property "power_in";
    interval ${RECORDINTERVAL};
	file "recorder_output.csv";
}
    ''')

    f.close()


def write_ESC_glm():
    # Write ESCs.glm file with the DERs/controllers selected in the config.json file (from frontend in GRIP platform)
    with open('config.json','r') as f:
        data = json.load(f)

    ESCs = data['controllers']

    f = open('ESCs.glm',"w+", newline='')
    for i in range(len(ESCs)):

        f.write("object meter {\n")
        f.write("    name ES_{}_meter;\n".format(i+1))
        f.write("    measured_energy_delta_timestep ${RECORDINTERVAL};\n")
        f.write("    phases {};\n".format(ESCs[i]['controller_phase']))
        f.write("    nominal_voltage 2401.7771;\n")
        f.write("}\n\n")

        f.write("object switch {\n")
        f.write("    name ES_{}_switch;\n".format(i+1))
        f.write("    from {};\n".format(ESCs[i]['node_name']))
        f.write("    to ES_{}_meter;\n".format(i+1))
        f.write("    phases {};\n".format(ESCs[i]['controller_phase']))
        f.write("    status CLOSED;\n")
        f.write("}\n\n")

        f.write("object inverter {\n")
        f.write("    name ES_{}_inv;\n".format(i+1))
        f.write("    phases {};\n".format(ESCs[i]['controller_phase']))
        f.write("    parent ES_{}_meter;\n".format(i+1))
        f.write("    rated_power {} kVA;\n".format(ESCs[i]['apparent_power_capacity_limit']))
        f.write("    generator_status ONLINE;\n")
        f.write("    generator_mode CONSTANT_PQ;\n")
        f.write("    inverter_efficiency 0.95;\n")
        f.write("    inverter_type FOUR_QUADRANT;\n")
        f.write("    four_quadrant_control_mode CONSTANT_PQ;\n")
        f.write("    P_Out 0.000 MW;\n")
        f.write("    Q_Out 0.000 MVAr;\n")
        f.write("}\n\n")

        f.write("object solar {\n")
        f.write("    name ES_{}_solar;\n".format(i+1))
        f.write("    generator_mode SUPPLY_DRIVEN;\n")
        f.write("    generator_status ONLINE;\n")
        f.write("    panel_type SINGLE_CRYSTAL_SILICON;\n")
        f.write("    parent ES_{}_inv;\n".format(i+1))
        f.write("    rated_power {} kW;\n".format(ESCs[i]['real_power_capacity_limit']))
        f.write("}\n\n")
    f.close()
