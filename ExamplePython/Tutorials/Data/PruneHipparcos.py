#!/usr/bin/python
import sys




def main():
  #removes all lines containing a non-numeric value
  f=open(sys.argv[1])
  o=open(sys.argv[1]+'.prune','w')
  for i,line in enumerate(f.readlines()):
    if i==0:
      o.write(line)
    else:
      outline=[]
      for item in line.split(','):
        if item.strip():
          try:
            float(item)
            outline+=[item]
          except:
            outline+=['NaN']
      o.write(','.join(outline))

if __name__ == '__main__':
  main()