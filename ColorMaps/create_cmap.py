#!/usr/bin/env python
class Color:
  def __init__(self,tup):
    self.r,self.g,self.b=[float(x)/255. for x in tup.split(',')]
    
if __name__ == '__main__':
  import sys
  name=sys.argv[1]
  colors=[Color(x) for x in sys.argv[2:]]
  print """<ColorMap name="%s" space="RGB">""" % name
  step_size=1/(len(colors)-1.)
  step=0
  for c in colors:
    print """<Point x="%f" o="1" r="%f" g="%f" b="%f"/>""" % \
      (step,c.r,c.g,c.b)
    step+=step_size
  print """</ColorMap>"""