import os

def replace_text(file_in,file_out,old_text,new_text):
    os.system("sed" + " 's|"+old_text+"|"+new_text+"|g'" + " " + file_in + " > " + file_out)

Compiler_old = "FC = ifx"
Compiler_new = "FC = ifx"

home=os.getcwd()

# install tools/
os.chdir("tools/FormCalc_Tools")
print("install FormCalc_Tools ...")
replace_text("./makefile","./makefile.tmp",Compiler_old,Compiler_new)
os.system("make -f makefile.tmp clean")
os.system("make -f makefile.tmp")
os.system("rm -f makefile.tmp")
os.chdir("../integrators/basesv5.1")
print("install basesv5.1 ...")
replace_text("./makefile","./makefile.tmp",Compiler_old,Compiler_new)
os.system("make -f makefile.tmp clean")
os.system("make -f makefile.tmp")
os.system("rm -f makefile.tmp")
os.system("rm -f *.o")
os.chdir("../../loops/LoopInts")
print("install LoopInts ...")
replace_text("./makefile","./makefile.tmp",Compiler_old,Compiler_new)
os.system("make -f makefile.tmp")
os.system("rm -f makefile.tmp")
os.chdir("../../utilities")
replace_text("./makefile","./makefile.tmp",Compiler_old,Compiler_new)
os.system("make -f makefile.tmp clean")
os.system("make -f makefile.tmp")
os.system("rm -f makefile.tmp")
os.system("rm -f *.o")
