from expdesign import *

dtype = 1

options = default.copy()
options['--experiment'] = 'mossbauer'
options['--type'] = '0'
options['--dim'] = '3'
options['--design'] = '-1.0,0,1.0' if dtype == 1 else '-2,0,2'
options['--sigeps'] = '0.1'

suffix = '';

options['--mcmcfile'] = 'expdesign/out/mcmc/design%d%s' % (dtype, suffix)
options['--nmcmc'] = 1000000

arguments = [executable]
for k, v in options.iteritems():
  arguments.append('%s' % k)
  arguments.append('%s' % v)


print ' '.join(arguments)

output = subprocess.check_output(arguments)
print output


