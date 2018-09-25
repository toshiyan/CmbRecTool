
import os

#//// text adds to a file ////#
def add(text,f,ini=False,opt=''):
  if ini:     os.system('echo '+opt+' "'+text+'" > '+f)
  if not ini: os.system('echo '+opt+' "'+text+'" >> '+f)

