import subprocess

# program_list = ["python3 andy.py -tao_monitor -tao_max_it 100 -l -0.40 -m '1_to_3_mesh.msh' -o '1_to_3_mesh_1' -er 0.1 -es 1.0",
#                 "python3 andy.py -tao_monitor -tao_max_it 100 -l -0.55 -m '1_to_3_mesh.msh' -o '1_to_3_mesh_2' -er 0.1 -es 1.0",
#                 "python3 andy.py -tao_monitor -tao_max_it 100 -l -0.60 -m '1_to_3_mesh.msh' -o '1_to_3_mesh_3' -er 0.1 -es 1.0"]

program_list = ["python3 andy.py -tao_monitor -tao_max_it 1000 -l 0.05 -m '1_to_3_mesh.msh' -o 'output1' -er 1.0 -es 1.0e-1 -k 8.0e-5 -e 4.0e-3 -p 1.0"]

i = 1
for program in program_list:
    print("-------------------------------------------------------")
    print("")
    print("Running test #{}".format(i))
    print("")
    print(program)
    print("")
    subprocess.run(program, shell = True)
    i = i + 1
