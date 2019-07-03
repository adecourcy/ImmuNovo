import os

if __name__ == "__main__":

  currentDirectory = os.getcwd()
  if " " in currentDirectory:
    print("\nUnfortunately, a dependency of this program will not compile")
    print("if there are any spaces in the directory path. This problem")
    print("can be fixed by moving this directory to a path without spaces,")
    print("running this script, and then moving the directory back\n")
    exit()

  os.system("unzip gc-7.6.0.zip")
  os.chdir("{}/gc-7.6.0".format(currentDirectory))
  os.system("./configure --prefix={} --disable-threads".format(currentDirectory))
  os.system("make")
  os.system("make check")
  os.system("make install")
  os.chdir(currentDirectory)
  os.system("cc -I{0}/include *.c {0}/lib/libgc.a -o findPep -lm -std=c99".format(currentDirectory))

