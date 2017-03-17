import os
import argparse
import sys

parser = argparse.ArgumentParser()
# Make the input file an input to the script
parser.add_argument("-i", "--input", type=str, help="The input file to be run over")

args = parser.parse_args()

print args.input

# Print errors if input file not given or not found
if not (args.input):
    print "ERROR"
    print "Needs input file"
    parser.print_help()
    sys.exit(1)

if not (os.path.isfile(args.input)):
    print "ERROR"
    print "input file %s not found" %(args.input)

# File for each efficiency (signal / background)
Ratios_File_Name = "Acceptance.dat"

# The peak values and left and right widths
SinglePeak = [(2599, 10, 10), (847, 10, 10), (1238, 10, 10), (511, 10, 10), (1771, 10, 10), (1037, 10, 10), (1360, 10, 10), (3202, 10, 10), (3451, 10, 10)]
DoublePeak = [(3254, (1.876 / 7.923), (3273.079 / 3253.503), 10, 35), (2035, (3.016 / 7.77), (2015.215 / 2034.791), 35, 10)]

# the filename for the tmp file to be created for sed
tmpfile = "cpp_files/tmpfile.cc"
# Final file for all the efficiencies
Efficiency_File = "Output/Efficiency.dat"
os.system("rm %s" %(Efficiency_File))

if not os.path.isdir("cpp_files"): os.system("mkdir cpp_files")
if not os.path.isdir("Output"): os.system("mkdir Output")

# Loop over double peaks
for (Name, Amplitude_Ratio, Mean_Ratio, Left_Window, Right_Window) in DoublePeak:

    Efficiency_File_Peak = "Output/Ratio_"+ str(Name) + ".dat"

    Cpp_Name = "cpp_files/peak_%s.cc" %(Name)
    
    os.system("cp doublepeak_template.cc %s" %(Cpp_Name))
# replace values with sed
    os.system("sed 's/_peak_number_/%s/g' <%s >%s; mv %s %s" %(Name, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
    os.system("sed 's/_amplitude_ratio_/%s/g' <%s >%s; mv %s %s" %(Amplitude_Ratio, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
    os.system("sed 's/_mean_ratio_/%s/g' <%s >%s; mv %s %s" %(Mean_Ratio, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
    os.system("sed 's/_left_window_/%s/g' <%s >%s; mv %s %s" %(Left_Window, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
    os.system("sed 's/_right_window_/%s/g' <%s >%s; mv %s %s" %(Right_Window, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
    os.system("sed 's;_filename_;%s;g' <%s >%s; mv %s %s" %(args.input, Cpp_Name, tmpfile, tmpfile, Cpp_Name))

    # run root and copy output to file
    os.system("root %s -l -q -b" %(Cpp_Name))
    os.system("cat %s >> %s" %(Efficiency_File_Peak, Efficiency_File))

# Do the same thing for single peaks
for (Name, Left_Window, Right_Window) in SinglePeak:
    
    Efficiency_File_Peak = "Output/Ratio_"+ str(Name) + ".dat"

    Cpp_Name = "cpp_files/peak_%s.cc" %(Name)
    
    os.system("cp singlepeak_template.cc %s" %(Cpp_Name))
    
    os.system("sed 's/_peak_number_/%s/g' <%s >%s; mv %s %s" %(Name, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
    os.system("sed 's/_left_window_/%s/g' <%s >%s; mv %s %s" %(Left_Window, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
    os.system("sed 's/_right_window_/%s/g' <%s >%s; mv %s %s" %(Right_Window, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
    os.system("sed 's;_filename_;%s;g' <%s >%s; mv %s %s" %(args.input, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
        
    os.system("root %s -l -q -b" %(Cpp_Name))

    os.system("cat %s >> %s" %(Efficiency_File_Peak, Efficiency_File))

# kill the tmpfile since we no longer need it
os.system("rm %s -f" %(tmpfile))

