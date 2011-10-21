#!/usr/bin/env python
# ----------------------------------------------------------------
# Gamma-Scout Decoder
#
# Version 0.1 - First released version (Lionel Bergeret)
# ----------------------------------------------------------------
# The contents of this file are distributed under the CC0 license.
# See http://creativecommons.org/publicdomain/zero/1.0/
# ----------------------------------------------------------------

import sys
import os
import struct
import datetime
from optparse import OptionParser

try:
    from pygooglechart import Chart
    from pygooglechart import SimpleLineChart
    from pygooglechart import Axis
    gChart = 1
except:
    gChart = 0
    
try:
    import matplotlib.pyplot as plt
    from matplotlib.dates import date2num
    from pylab import *
    from numpy import *
    pChart = 1
except:
    pChart = 0

def Hex2Float(h):
    # Dr.Mirow formula
    # z = m when e = 0
    # z = (2^10 + m)*2^(e-1) when e > 0

    v = struct.unpack('>H', h.decode('hex'))[0]

    e = int((v >> 10) & 0x0000001f)    # 6 bits exponent
    m = int(v & 0x000003ff)            # 10 bits mantissa
 
    e = (e+1)/2 # the two first bit are concatenate to one (Dr.Mirow)

    if (e == 0):
      return m
    else:
      return (1024 + m)*2**(e-1)
      
def EMA(s, n):
    """
    returns an n period exponential moving average for
    the time series s

    s is a list ordered from oldest (index 0) to most recent (index
    -1) n is an integer

    returns a numeric array of the exponential moving average
    """
    s = array(s)
    ema = []
    j = 1
    #get n sma first and calculate the next n period ema
    sma = sum(s[:n]) / n
    multiplier = 2 / float(1 + n)
    ema.append(sma)
    #EMA(current) = ( (Price(current) - EMA(prev) ) xMultiplier) + EMA(prev)
    ema.append(( (s[n] - sma) * multiplier) + sma)
    #now calculate the rest of the values
    for i in s[n+1:]:
        tmp = ( (i - ema[j]) * multiplier) + ema[j]
        j = j + 1
        ema.append(tmp)
    return ema

def SMA(s, n):
    """
    returns an n period moving average for the time series s

    s is a list ordered from oldest (index 0) to most recent (index -1)
    n is an integer

        returns a numeric array of the moving average

    See also ema in this module for the exponential moving average.
    """
    s = array(s)
    c = cumsum(s)
    return (c[n-1:] - c[:-n+1]) / float(n-1)
    
def FillOver(ax, x, y, val, color, over=True):
    """
    Plot filled x,y for all y over val
    if over = False, fill all areas < val
    """
    ybase = asarray(y)-val
    crossings = nonzero(less(ybase[:-1] * ybase[1:],0))

    if ybase[0]>=0: fillon = over
    else:           fillon = not over


    indLast = 0
    for ind in crossings:
        if fillon:
            thisX = x[indLast:ind+1]
            thisY = y[indLast:ind+1]
            thisY[0] = val
            thisY[-1] = val
            ax.fill(thisX, thisY, facecolor=color)
        fillon = not fillon
        indLast = ind
      
def GammascoutDecode(filename):
    """
    returns a decoded dictionary from RAW data in filename
    """
    # Initialize raw data list
    gslog = []
    
    # Initialize decoded data dictionary
    decoded = {}
    decoded["time"]=[]
    decoded["count"]=[]
    decoded["interval"]=[]
    decoded["cps"]=[]
    decoded["cpm"]=[]
    decoded["error"]=[]
    
    # Open the raw log file
    f = open(filename, 'r')

    # Process the file
    while True:
       line = f.readline()
   
       if len(line) == 0:
            break # EOF

       # Cleanup the log line
       line = line.lstrip(" ").rstrip()

       # process the valid lines
       if (len(line)>1):
         items = line.split(" ")
         if len(items) > 16:
           # remove the address
           items = items[1:] 
     
           for data in items:
             gslog.append(data)

    # Sample record
    # 0100 fe 10 21 28 07 11 f1 ff 01 00 00 0f f2 03 5e 03
    # 0110 7e 03 b0 03 59 03 5c 03 9b 03 89 03 6f 03 85 03
    # 0120 64 03 b0 03 54 03 9b 03 dc 03 75 03 8f 03 97 03
    # 0130 a1 03 75 03 a5 03 5b 03 50 03 6a 03 6a 03 7f 03
    # 0140 6b 03 74 03 7d 03 7a 03 97 03 88 03 a3 03 25 03
    # 0150 57 03 5c 03 5c 03 27 03 8e 03 29 03 55 03 4d 0d
    # 0160 43 03 6f 03 99 03 8b 02 d4 02 66 02 1e 02 20 02
    # 0170 3b 02 73 06 3f 03 9a 03 7a 03 a2 03 91 03 63 03
    # 0180 8f 03 7a 03 7a 03 91 03 78 03 ad 03 76 03 96 03
    # 0190 96 05 41 03 8d 03 79 03 98 03 ad 03 7e 03 96 03
    # 01a0 5e 03 68 03 28 03 41 03 20 03 3f 03 2f 03 45 03
    # 01b0 20 02 f4 02 db 02 30 02 3e 02 32 02 50 02 48 02
    # 01c0 61 02 65 02 7c 02 b4 03 02 03 76 03 51 03 67 03
    # 01d0 55 03 52 03 4d 03 1a 03 51 03 43 03 6e 03 69 03
    # 01e0 3c 03 4c 02 ab 02 be 02 63 02 4d 02 8d 02 60 02
    # 01f0 77 02 82 02 75 02 c9 ff ff ff ff ff ff ff ff ff
    
    # TODO: Check if it is really correct
    timeDictionary = {"f0": 7*24*60*60, "f1": 24*60*60 , "f2": 60*60, "f3": 10*60, "f4": 60}

    print "Gamma Scout ID=%s%s%s" % (gslog[2], gslog[1], gslog[0])

    # Initialize timestamp
    timestamp = datetime.datetime.now()
    
    # Initialize intervals
    previousIntervalSec = 0
    newIntervalSec = 0
    timeDelta = datetime.timedelta(0)

    # Retrieve the data end address
    finalAddress = Hex2Float(gslog[33]+gslog[32])

    # Decode the binary data
    i = 256 # go to address 0x100 directly
    while i < finalAddress:
      # ----------------------------------------------------------------
      # Start time definition
      # ----------------------------------------------------------------
      if (gslog[i] == "fe") and (gslog[i+6] in timeDictionary.keys()):
        # New start time
        year = int("20"+gslog[i+5])
        month = int(gslog[i+4])
        day = int(gslog[i+3])
        hour = int(gslog[i+2])
        minute = int(gslog[i+1])
        
        # Interval definition
        newIntervalSec = timeDictionary[gslog[i+6]]
        timeDelta = datetime.timedelta(seconds=newIntervalSec)

        # Timestamp definition
        timestamp = datetime.datetime(year, month, day, hour, minute, 0)

        i=i+7 # skip the 7 processed bytes
        
      # ----------------------------------------------------------------
      # Interval transition with counts
      # TODO: Check if it is correct
      # ----------------------------------------------------------------
      elif (gslog[i] == "ff") and (i < len(gslog)-5) and (gslog[i+5] in timeDictionary.keys()):
        previousIntervalSec = int(gslog[i+1], 16) * 60 + int(gslog[i+2], 16) * 60 *60
        previousCounts = int(gslog[i+3]+gslog[i+4], 16)
        newIntervalSec = timeDictionary[gslog[i+5]]
        if previousIntervalSec == 0:
          previousIntervalSec = newIntervalSec
      
        timeDelta = datetime.timedelta(seconds=previousIntervalSec)
        cps = float(previousCounts)/float(previousIntervalSec)
        cpm = float(previousCounts)/float(previousIntervalSec/60)

        if previousCounts > 0:  # Ignore 0 measures
          decoded["time"].append("%s" % timestamp)
          decoded["count"].append(previousCounts)
          decoded["interval"].append(previousIntervalSec)
          decoded["cps"].append(cps)
          decoded["cpm"].append(cpm)
          decoded["error"].append(100.0/sqrt(previousCounts))

        timestamp += timeDelta
    
        # New Interval definition
        timeDelta = datetime.timedelta(seconds=newIntervalSec)
        i=i+6 # skip the 6 processed bytes 
      
      # ----------------------------------------------------------------
      # New Interval definition
      # ----------------------------------------------------------------
      elif (gslog[i] in timeDictionary.keys()):
        newIntervalSec = timeDictionary[gslog[i]]
        timeDelta = datetime.timedelta(seconds=newIntervalSec)
        i=i+1 # skip 1 byte
        
      # ----------------------------------------------------------------
      # Measurements
      # ----------------------------------------------------------------
      elif (gslog[i] != "ff") and (gslog[i] != "20"):
        data = Hex2Float(gslog[i]+gslog[i+1])
        cps = float(data)/float(newIntervalSec)
        cpm = float(data)/float(newIntervalSec/60)
        
        if data > 0: # Ignore 0 measures
          decoded["time"].append("%s" % timestamp)
          decoded["count"].append(data)
          decoded["interval"].append(newIntervalSec)
          decoded["cps"].append(cps)
          decoded["cpm"].append(cpm)
          decoded["error"].append(100.0/sqrt(data))

        timestamp += timeDelta
        i=i+2 # skip the 2 processed bytes 
        
      # ----------------------------------------------------------------
      # Unknown byte
      # ----------------------------------------------------------------
      else:
        i=i+1 # skip 1 byte
  
    f.close()
    return decoded
    
def RenderCSV(data, filename):
    """
    create a CSV file from decoded data under filename (.csv)
    """
    print "Rendering CSV [%s]" % filename
    f = open(filename, "w")
    f.write("date, interval, counts, cps, cpm, error\n");
    for i in range(len(data["time"])):
        f.write("%s, %d, %d, %f, %f, %f\n" % (data["time"][i], data["interval"][i], data["count"][i], data["cps"][i], data["cpm"][i], data["error"][i]))
    f.close
    
def RenderGnuplot(data, filename):
    """
    create a GnuPlot space separated data file from decoded data under filename (.data)
    """
    print "Rendering Gnuplot [%s]" % filename
    f = open(filename, "w")
    for i in range(len(data["time"])):
        f.write("%s %d %f %f\n" % (data["time"][i], data["count"][i], data["cps"][i], data["cpm"][i]))
    f.close
        
def RenderGoogleChart(data, filename):    
    """
    create a GoogleChart from decoded data under filename (.png)
    """
    print "Rendering GoogleChart [%s]" % filename
    # Retrieve chart data
    elements=[]
    max_y = 0
    min_y = 9999
    for i in range(len(data["time"])):
        if data["cps"][i] > max_y: max_y = data["cps"][i]
        if data["cps"][i] < min_y: min_y = data["cps"][i]
        elements.append(data["cps"][i])

    # Chart size of 600x375 pixels and specifying the range for the Y axis
    chart = SimpleLineChart(600, 375, y_range=[min_y-0.1, max_y+0.1])
    
    # Add the chart data
    chart.add_data(elements)
    
    # Set the line colour to blue
    chart.set_colours(['0000FF'])

    # Set the vertical stripes
    chart.fill_linear_stripes(Chart.CHART, 0, 'CCCCCC', 0.2, 'FFFFFF', 0.2)

    # Set the horizontal dotted lines
    chart.set_grid(0, 25, 5, 5)

    # Define the Y axis labels
    left_axis = [x * 0.1 for x in range(0, int(max_y/0.1))]
    left_axis[0] = 'CPS'
    chart.set_axis_labels(Axis.LEFT, left_axis)

    chart.download(filename)

hlcount = 0
def Highlight(time, data, index, unit, reason):
    global hlcount
    if index and index > 0:
      x,y = (time[index], data[index])
      plot([x], [y], 'ko')
      annotate('%0.3f %s\n%s\n(%s)' % (y, unit, num2date(time[index]).strftime('%Y-%m-%d %H:%M'), reason), xy=(x,y),
         xytext=(20, -(hlcount+4.7)*22), textcoords='axes points',
         arrowprops=dict(arrowstyle="-",
         connectionstyle="angle,angleA=0,angleB=80,rad=10"),
         horizontalalignment='left',
         verticalalignment='bottom',
         fontsize=8)
      hlcount+=1
    
def RenderMatplotlib(data, filename):
    """
    create a Matplotlib chart from decoded data under filename (.pdf)
    """
    print "Rendering Matplotlib [%s]" % filename
    
    # Create the figure (A4 format)
    plt.figure(num=None, figsize=(8.27, 11.69), dpi=100)
    unixtime = date2num( [datetime.datetime.strptime(x, "%Y-%m-%d %H:%M:%S") for x in data['time']] )

    # ------------------------------------------------------------------
    # First graph uSievert
    # ------------------------------------------------------------------
    ax1 = plt.subplot2grid((4, 2), (0, 0), rowspan=2, colspan=2)
    microSievert = [x/123.0 for x in data['cpm']]

    plt.plot(unixtime, microSievert, 'b-', label='uSv/h') # alpha=0.4

    # Compute the trend line
    z = polyfit(unixtime, microSievert, 1)
    p = poly1d(z)
    plt.plot(unixtime, p(unixtime), "r--", label='Trend')
    
    sma = SMA(microSievert, 8)
    plt.plot(unixtime[8-1:], sma, 'g--', linewidth=1, label='SMA 8pt')
    
    #ema = EMA(microSievert, 8)
    #plt.plot(unixtime[8-1:], ema, 'y--', linewidth=1, label='EMA 8pt')
    #FillOver(plt, unixtime, microSievert,  0.2,  "purple", over=True)
    
    # Highlight maximum
    hlcount = 0
    index = np.argmax(microSievert)
    Highlight(unixtime, microSievert, index, "uSv/h", "maximum")
    
    plt.legend(loc='upper left', fancybox=True, shadow=True, prop=dict(size=8))
    plt.grid(True, which="both", linestyle="dotted")
    plt.xlabel("Date & Time", fontsize=7)
    plt.ylabel("microSieverts Per Hour", fontsize=7)
    plt.xticks(fontsize=7, rotation=90)
    plt.yticks(fontsize=7)

    # ------------------------------------------------------------------
    # Second graph CPM
    # ------------------------------------------------------------------
    ax2 = plt.subplot2grid((4, 2), (2, 0), colspan=2, sharex=ax1)
    plt.plot(unixtime, data['cpm'], 'b-', label='CPM')
    
    # Compute the trend line
    z = polyfit(unixtime, data['cpm'], 1)
    p = poly1d(z)
    plt.plot(unixtime, p(unixtime), "r--", label='Trend')
    
    plt.legend(loc='upper left', fancybox=True, shadow=True, prop=dict(size=8))
    plt.grid(True, which="both", linestyle="dotted")
    plt.xlabel("Date & Time", fontsize=7)
    plt.ylabel("Counts Per Minute", fontsize=7)
    plt.xticks(fontsize=7, rotation=90)
    plt.yticks(fontsize=7)
    
    # ------------------------------------------------------------------
    # Third graph statistical errors
    # ------------------------------------------------------------------
    ax3 = plt.subplot2grid((4, 2), (3, 0),  colspan=2, sharex=ax1)
    plt.errorbar(unixtime, data['cpm'], yerr=data['error'], label='Percent', marker='.', linestyle='-', fmt='ro')
    
    plt.legend(loc='upper left', fancybox=True, shadow=True, prop=dict(size=8))
    plt.grid(True, which="both", linestyle="dotted")
    plt.xlabel("Date & Time", fontsize=7)
    plt.ylabel("Statistical error", fontsize=7)
    plt.xticks(fontsize=7, rotation=90)
    plt.yticks(fontsize=7)

    # ------------------------------------------------------------------
    # Save to PDF
    # ------------------------------------------------------------------
    # Set major x ticks
    # TODO: add command line option !!!
    #ax1.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=6)) # Every 6 hours
    #ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m/%d\n%H:%M'))
    ax1.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=24)) # Daily
    ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m/%d'))

    plt.savefig(filename)
    #plt.show()
    
# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
if __name__=="__main__":
    parser = OptionParser("Usage: GammaScoutDecode [options] <raw-data-file>")
    parser.add_option("-f", "--file", dest="filename",
                      help="write report to FILE", metavar="FILE")
    parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="don't print status messages to stdout")

    (options, args) = parser.parse_args()
    
    if len(args) != 1:
        parser.error("wrong number of arguments")
        
    decoded = GammascoutDecode(args[0])
    RenderCSV(decoded, args[0]+".csv")
    RenderGnuplot(decoded, args[0]+".data")
    #if gChart:
    #  RenderGoogleChart(decoded, args[0]+".png")
    if pChart:
      RenderMatplotlib(decoded, args[0]+".pdf")

