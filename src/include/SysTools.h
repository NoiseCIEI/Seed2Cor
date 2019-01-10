#ifndef SYSTOOLS_H
#define SYSTOOLS_H
#include <vector>
#include <string>

/* prompt and continue in 'time' seconds if no input is received */
void TimedContinue (int time);

/* read and compute the currently available RAM from file '/proc/meminfo' */
void EstimateMemAvail (long &MemAvail);

/* create a directory named '*dirname', exit upon failure
 * returns true if dir is made successfully and false if already exist */
bool MKDir(const char *dirname);
bool MKDirs(const char *dirname);

/* delete the file (or empty dir) named '*fname', prompt to continue upon failure */
void fRemove (const char *fname);

/* descend into the directory '*dirname' and remove everything recursively */
int dRemove(const char *dirname);

/* rename a file or a directory, prompt to continue upon failure */
void Move (const char *oldname, const char *newname);

/* move any file in '*odir' that matches '*pattern' into '*tdir and set '*nfile' to the number 
   of files moved return a list of moved file names if retlst==1 and return NULL otherwise */ 
//char * wMove (const char *odir, const char *pattern, const char *tdir, int retlst, int *nfile);
bool wMove (const char *odir, const char *pattern, const char *tdir, std::vector<std::string> &outlist);

/* make a copy of file '*oldname' to '*newname', prompt to continue upon failure */
void Copy(const char *oldname, const char *newname);

/* return a list of file names from '*dir' that matche '*pattern'
   the input 'type' parameter decides how sub-directories are handdled:
   0: list files int the root dir only
   1: list files and directories in the root dir
   2: list all files
   3: list all files and directories */
//char * List(const char *dir, const char *pattern, int type, int *nfile);
bool List(const char *dir, const char *pattern, int type, std::vector<std::string> &filelist);

#endif
