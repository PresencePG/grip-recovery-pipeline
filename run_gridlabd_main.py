import gridlabd
import traceback
import os
import sys, getopt


path = ''
modelinputfile = ''
outputfile = ''
compile_flag=''



def command_line_template():
    sys.stderr.write('run_gridlabd_main.py -C -W <inputfilepath> -i <modelinputfile> -o <outputfile>')
    sys.stderr.write('-C : [OPTIONAL] Compile flag for GridLAB-D.')
    sys.stderr.write('-W : [OPTIONAL] Specifies folder location of GLM file, if not specified assumes current folder.')
    sys.stderr.write('-i : [REQUIRED] Name of the model input file')
    sys.stderr.write('-o : [OPTIONAL] JSON dump file name.')
    return True


try:
    try:
        opts, args = getopt.getopt(sys.argv[1:], "CW:i:o:", ["Wfile=", "ifile=", "ofile="])
    except getopt.GetoptError:
        sys.stderr.write('Incorrect command arguments')
        command_line_template()
        sys.exit(2)
    if not opts:
        sys.stderr.write('Missing command arguments')
        command_line_template()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-C':
            compile_flag = 1
        elif opt in ("-W", "--Wfile"):
            path = arg
        elif opt in ("-i", "--ifile"):
            modelinputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        else:
            print("ERROR: Incorrect argument set: python3 [file] -W [path] -i [modelinputfile] -o [outputfile]")
    if path:
        print("Setting path to folder: ", path)
        os.chdir(path)
    else:
        print("Setting path to current folder")

    gridlabd.command('--profile')
    gridlabd.command('--version')
    gridlabd.command('--relax')
    if compile_flag == 1:
        print("Compile Flag is Set")
        gridlabd.command('-C')
    else:
        import GLM_writer
        # Create updated recovery_ieee37.GLM
        # GLM_writer.write_recovery_ieee37(modelinputfile)
        GLM_writer.write_glm_flex(modelinputfile)
        # Create ESCs.GLM
        GLM_writer.write_ESC_glm()

    gridlabd.command('--redirect')
    gridlabd.command('all')
    gridlabd.command(modelinputfile)
    if outputfile:
        gridlabd.command('-o')
        gridlabd.command(outputfile)
    gridlabd.start('wait')
    print('\n')
    if os.stat('gridlabd.err').st_size == 0:
        print('No errors in file.')
    else:
        sys.stderr.write('Error(s) found in the model, printing gridlabd.err file contents...')
        fr_err = open('gridlabd.err', 'r')
        sys.stderr.write(fr_err.read())
        sys.exit(3)
    if os.stat('gridlabd.wrn').st_size == 0:
        print('No errors in file.')
    else:
        print('Warning(s) found in the model, printing gridlabd.wrn file contents...')
        fr_wrn = open('gridlabd.wrn', 'r')
        print(fr_wrn.read())


except Exception as e:
    sys.stderr.write(f'Unhandled Exception in main.py: {e}')
    sys.stderr.write(traceback.format_exc())
    raise e
