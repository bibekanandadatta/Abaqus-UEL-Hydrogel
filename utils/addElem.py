input_fname = input("ABAQUS input file: ")
elem        = input("Enter element type (1. TET4, 2. HEX8, 3. TRI3-AX, 4. QUAD4-AX, 5. TRI3-PE, 6. QUAD4-PE): ")
jtype       = int(elem)


offset      = 100000

if jtype == 1:
    oldElemStr  = '*Element, type=C3D4'
    newElemStr  = '*User Element,Type=U1,Nodes=<nnode>,Coordinates=<nDim>,Properties=<matProps>,Iproperties=<intProps>,Variables=<nsvars>,Unsymm\n'\
                  '1,2,3,11\n'\
                  '*Element, type=U1'
elif jtype == 2:
    oldElemStr  = '*Element, type=C3D8'
    newElemStr  = '*User Element,Type=U2,Nodes=<nnode>,Coordinates=<nDim>,Properties=<matProps>,Iproperties=<intProps>,Variables=<nsvars>,Unsymm\n'\
                  '1,2,3,11\n'\
                  '*Element, type=U2'

elif jtype == 3:
    oldElemStr  = '*Element, type=CAX3'
    newElemStr  = '*User Element,Type=U3,Nodes=<nnode>,Coordinates=<nDim>,Properties=<matProps>,Iproperties=<intProps>,Variables=<nsvars>,Unsymm\n'\
                  '1,2,11\n'\
                  '*Element, type=U3'
elif jtype == 4:
    oldElemStr  = '*Element, type=CAX4'
    newElemStr  = '*User Element,Type=U4,Nodes=<nnode>,Coordinates=<nDim>,Properties=<matProps>,Iproperties=<intProps>,Variables=<nsvars>,Unsymm\n'\
                  '1,2,11\n'\
                  '*Element, type=U4'
elif jtype == 2:
    oldElemStr  = '*Element, type=CPE3'
    newElemStr  = '*User Element,Type=U3,Nodes=<nnode>,Coordinates=<nDim>,Properties=<matProps>,Iproperties=<intProps>,Variables=<nsvars>,Unsymm\n'\
                  '1,2,11\n'\
                  '*Element, type=U3'
elif jtype == 4:
    oldElemStr  = '*Element, type=CPE4'
    newElemStr  = '*User Element,Type=U4,Nodes=<nnode>,Coordinates=<nDim>,Properties=<matProps>,Iproperties=<intProps>,Variables=<nsvars>,Unsymm\n'\
                  '1,2,11\n'\
                  '*Element, type=U4'
else:
    print("element is not supported")

# extra element to activate COUPLED TEMP-DISPLACEMENT step
extraElemStr    = '\n** EXTRA ELEMENT\n'\
                  '*Node\n'\
                  '999996, 0.0, 0.0\n'\
                  '999997, .00001, 0.0\n'\
                  '999998, .00001, .00001\n'\
                  '999999, 0.0, .00001\n'\
                  '*Element, Type=CPE4T\n'\
                  '999999,999996,999997,999998,999999\n'\
                  '*Nset, nset=extraElement\n'\
                  '999996,999997,999998,999999\n'\
                  '*Elset,elset=extraElement\n'\
                  '999999'

elem_info   = list()
elem_record = False


#replacing the original element type definition with new definition
with open(input_fname, "r") as fin:
    new_file_content = ""
    for line in fin:
        stripped_line = line.strip()
        new_line = stripped_line.replace(oldElemStr,newElemStr)
        new_file_content += new_line +"\n"
    fin.close()

with open(input_fname, "w") as fout:
    fout.write(new_file_content)
    fout.close()

#reading original element connectivity for dummy elements
with open(input_fname, 'r') as fin:
    lines = fin.readlines()
    for line in lines:
        line_splt = line.strip().split(',')
        
        if len(line_splt) and line_splt[0] == '*Element':
            elem_record = True
            continue
        
        if len(line_splt) and line_splt[0] == '*System':
            elem_record = False
        
        if elem_record is True:
            elem_list = [int(x) for x in line_splt]
            elem_list[0] += offset
            elem_info.append(elem_list)
    fin.close()

with open(input_fname,'a+') as fout:

    #add dummy elements to the file
    if jtype == 1:
        fout.write('*Element, type=C3D4')
    elif jtype == 2:
        fout.write('*Element, type=C3D8')
    elif jtype == 3:
        fout.write('*Element, type=CAX3')
    elif jtype == 4:
        fout.write('*Element, type=CAX4')
    elif jtype == 5:
        fout.write('*Element, type=CPE3')
    elif jtype == 6:
        fout.write('*Element, type=CPE4')
    else:
        print("element is not supported")

    for elem_list in elem_info:
        fout.write('\n')
        fout.write(str(elem_list)[1:-1])
    fout.write('\n*Elset, elset=ElDummy, generate\n')
    fout.write(str(elem_info[0][0]) + "," + str(elem_info[-1][0]) + "," + str(1))

    #add extra element block to the input file
    fout.write(extraElemStr)

    fout.close()