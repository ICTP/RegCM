* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* use.gs
* 
* This function is similar to the GrADS command 'open'. It checks 
* the list of current files to see if the specified file has already 
* been opened. If there is a match, then the default file is changed 
* to the matching file. If the specified file is not a current file,
* then it is opened and set to the default file.
*
* Written by J.M.Adams, May 2001, updated July 2001
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
function use(args)

* Only argument is the name of the file to open
  if (args!='') 
    openfile = subwrd(args,1)
  else 
    prompt '---> Enter the name of the file you want to use: '
    pull openfile
  endif

* Add ".ctl" to the filename if it's not already there
  len  = math_strlen(openfile)
  if (len >= 4)
    tail = substr(openfile,len-3,4)
  else 
    tail = openfile
  endif
  if (tail != '.ctl')
    openfile = openfile%'.ctl'
  endif

* Get a list of all open files 
  'q files'
  filelist = result

* If no files are open, 'openfile' becomes File 1
  if (subwrd(filelist,1) = 'No')
    'open 'openfile
    line1 = sublin(result,1)
    line2 = sublin(result,2)
    if (subwrd(line2,2) = 'Error:')
       say openfile' not opened'
    else
      'set dfile 1'
      prompt openfile' has been opened as File 1'
      say ' and is the default file'
    endif
  else 
*   Compare 'openfile' to all open files
    count = 0
    while(1)
      line1 = sublin(filelist,count*3+1)
      line2 = sublin(filelist,count*3+2)
      line4 = sublin(filelist,count*3+4)
      fnum = subwrd(line1,2)
      filename.fnum   = subwrd(line2,2)

*     If 'openfile' matches an existing filename, set a new default file
      if (filename.fnum = openfile)
        'set dfile 'fnum
        prompt openfile' is already open as File 'fnum
        say ' and is now the default file'
        break
      endif

*     If 'openfile' doesn't match any existing files, open it
      if line4 = ''
        'open 'openfile
        line1 = sublin(result,1)
        line2 = sublin(result,2)
        if (subwrd(line2,2) = 'Error:')
          say openfile' not opened'
        else
          newfnum = subwrd(line2,8)
          if (newfnum = fnum+1)
            'set dfile 'newfnum
            prompt openfile' has been opened as File '%(fnum+1)
            say ' and is the new default file'
          else
            say 'logic error on open for 'openfile
          endif
        endif
        break
      endif

      count = count + 1
    endwhile     
  endif

* THE END
