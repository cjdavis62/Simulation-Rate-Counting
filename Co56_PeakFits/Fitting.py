import os

Ratios_File_Name = "Acceptance.dat"


SinglePeak = [(2599, 10, 10), (847, 10, 10), (1238, 10, 10), (511, 10, 10), (1771, 10, 10), (1037, 10, 10), (1360, 10, 10), (3202, 10, 10), (3451, 10, 10)]
DoublePeak = [(3254, (1.876 / 7.923), (3273.079 / 3253.503), 10, 35), (2035, (3.016 / 7.77), (2015.215 / 2034.791), 35, 10)]


tmpfile = "cpp_files/tmpfile.cc"
Efficiency_File = "Output/Efficiency.dat"
os.system("rm %s" %(Efficiency_File))

if not os.path.isdir("cpp_files"): os.system("mkdir cpp_files")
if not os.path.isdir("Output"): os.system("mkdir Output")

for (Name, Amplitude_Ratio, Mean_Ratio, Left_Window, Right_Window) in DoublePeak:

    Efficiency_File_Peak = "Output/Ratio_"+ str(Name) + ".dat"

    Cpp_Name = "cpp_files/peak_%s.cc" %(Name)
    
    os.system("cp doublepeak_template.cc %s" %(Cpp_Name))

    os.system("sed 's/_peak_number_/%s/g' <%s >%s; mv %s %s" %(Name, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
    os.system("sed 's/_amplitude_ratio_/%s/g' <%s >%s; mv %s %s" %(Amplitude_Ratio, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
    os.system("sed 's/_mean_ratio_/%s/g' <%s >%s; mv %s %s" %(Mean_Ratio, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
    os.system("sed 's/_left_window_/%s/g' <%s >%s; mv %s %s" %(Left_Window, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
    os.system("sed 's/_right_window_/%s/g' <%s >%s; mv %s %s" %(Right_Window, Cpp_Name, tmpfile, tmpfile, Cpp_Name))

    os.system("root %s -l -q -b" %(Cpp_Name))
    os.system("cat %s >> %s" %(Efficiency_File_Peak, Efficiency_File))

for (Name, Left_Window, Right_Window) in SinglePeak:
    
    Efficiency_File_Peak = "Output/Ratio_"+ str(Name) + ".dat"

    Cpp_Name = "cpp_files/peak_%s.cc" %(Name)
    
    os.system("cp singlepeak_template.cc %s" %(Cpp_Name))
    
    os.system("sed 's/_peak_number_/%s/g' <%s >%s; mv %s %s" %(Name, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
    os.system("sed 's/_left_window_/%s/g' <%s >%s; mv %s %s" %(Left_Window, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
    os.system("sed 's/_right_window_/%s/g' <%s >%s; mv %s %s" %(Right_Window, Cpp_Name, tmpfile, tmpfile, Cpp_Name))
        
    os.system("root %s -l -q -b" %(Cpp_Name))

    os.system("cat %s >> %s" %(Efficiency_File_Peak, Efficiency_File))

os.system("rm %s -f" %(tmpfile))

