# ianh@astro.ox.ac.uk

import sys
import glob
import pickle
import numpy
import time
import datetime
from optparse import OptionParser
from pyrap.tables import table
from astropy.time import Time


def ms_info(msname):
	tt = table(msname,ack=False)
	uvw = tt.getcol('UVW')
	uv_max = numpy.max(((uvw[:,0]**2.0+uvw[:,1]**2.0)**0.5))
        times = tt.getcol('TIME')
        t_int = numpy.mean(tt.getcol('EXPOSURE'))
	tt.done()
	tt = table(msname+'/SPECTRAL_WINDOW',ack=False)
	freqs = tt.getcol('CHAN_FREQ')[0]
	tt.done()
	return uv_max,freqs,t_int,times


def get_flags(msname,uv0,uv1,t0,t1):
	# uv0,uv1 = uv range
	# t0,t1 = time range
	tt = table(msname,ack=False)
	taql = 'SQRT(SUMSQR(UVW[:2])) >= '+str(uv0)+' && SQRT(SUMSQR(UVW[:2])) < '+str(uv1)+' '
	taql += '&& TIME >= '+str(t0)+' && TIME <= '+str(t1)
	subtab = tt.query(query=taql)
	flags = subtab.getcol('FLAG')
	subtab.done()
	tt.done()
	return flags


parser = OptionParser(usage='%prog [options] msname \n Intervals of zero use entire range')
parser.add_option('-u','--uv',dest='delta_uv',help='Baseline averaging interval [m] (default=40)',default=40)
parser.add_option('-f','--freq',dest='delta_freq',help='Frequency averaging interval [chans] (default=4)',default=4)
parser.add_option('-t','--time',dest='delta_time',help='Time averaging interval [s] (default=0)',default=0)
parser.add_option('-o','--out',dest='opfile',help='Output file for flagging stats, created inside MS',default='flag_stats.dat')

(options,args)  = parser.parse_args()

delta_uv = options.delta_uv
delta_freq = options.delta_freq
delta_time = options.delta_time

if len(args) != 1:
	print 'Please specify a Measurement Set'
	sys.exit(-1)
else:
	msname = args[0].rstrip('/')

opfile = msname+'/'+options.opfile

now = str(datetime.datetime.now()).replace(' ','-').replace(':','-').split('.')[0]
clock0 = time.time()

print 'Measurement Set .............. '+msname

uv_max,freqs,t_int,times = ms_info(msname)

if delta_uv > 0:
	uv_bins = int(numpy.ceil(uv_max/delta_uv))
else:
	uv_bins = 1
if delta_freq > 0:
	freq_bins = int(numpy.ceil(len(freqs)/delta_freq))
else:
	freq_bins = 1
if delta_time > 0:
	time_bins = int(numpy.ceil(len(times)/delta_time))
else:
	time_bins = 1

#sys.exit(-1)

print ''
print 'Mean t_int ................... '+str(round(t_int,2))+' s'
print 'Number of time bins .......... '+str(time_bins)
print 'delta t ...................... '+str(delta_time)+' s'
print ''
print 'Max baseline ................. '+str(round(uv_max,2))+' m'
print 'Number of u,v bins ........... '+str(uv_bins)
print 'delta u,v .................... '+str(delta_uv)+' m'
print ''
print 'Number of chans .............. '+str(len(freqs))
print 'Number of chan bins .......... '+str(freq_bins)
print 'delta freq ................... '+str(delta_freq)
print ''

f = open(opfile,'w')
print >>f,'#generated='+now
print >>f,'#msname='+msname
print >>f,'#uv_bins='+str(uv_bins)
print >>f,'#freq_bins='+str(freq_bins)
print >>f,'#time_bins='
print >>f,'#delta_time='+str(delta_time)
print >>f,'#uv0,uv1,ch0,ch1,f0,f1,flag_frac,n_meas,t_mean,t_iso'
for t in range(0,time_bins):
	t0 = t*delta_time
	t1 = (t+1)*delta_time
	for i in range(0,uv_bins):
		print 'Processing chunk ............. '+str(i)+' of '+str(uv_bins)
		uv0 = i*delta_uv
		uv1 = (i+1)*delta_uv
		t0 = times[0]
		t1 = times[-1]
		t_mean = numpy.mean((t0,t1))
		t_obj = Time(t_mean/86400.0,format='mjd',scale='utc')
		t_iso = t_obj.iso.replace(' ','_')
		flags = get_flags(msname,uv0,uv1,t0,t1)
		for j in range(0,freq_bins):
			ch0 = (j*delta_freq)
			ch1 = ((j+1)*delta_freq) - 1
			f0 = freqs[ch0]
			f1 = freqs[ch1]
			if flags is not None:
				vals,counts = numpy.unique(flags[:,ch0:ch1,:],return_counts=True)
				n_meas = numpy.sum(counts)
				if len(vals) == 1 and vals == True:
					flag_frac = 1.0
				elif len(vals) == 1 and vals == False:
					flag_Frac = 0.0
				else:
					flag_frac = round(float(counts[1])/float(numpy.sum(counts)),3)
			else:
				flag_frac = 0.0
			print >>f,uv0,uv1,ch0,ch1,f0,f1,t0,t1,flag_frac,n_meas,t_mean,t_iso

print 'Wrote log file ............... '+opfile
elapsed = time.time() - clock0
print 'Elapsed time ................. '+str(round(elapsed/60.0,2))+' min'
print 'Done'
