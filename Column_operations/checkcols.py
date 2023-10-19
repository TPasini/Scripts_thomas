import casacore.tables
import os, lib_ms, lib_util, lib_log, glob

logger_obj = lib_log.Logger('testcol.log')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry = False)
MSs = lib_ms.AllMSs(glob.glob('TC00.MS'), s)

os.system('taql "ALTER TABLE TC00.MS DELETE COLUMN TEST"')
logger.info('Adding TEST column to file..')
MSs.run('addcol2ms.py -m TC00.MS -c TEST', log='testlog.log', commandType='python')

t = casacore.tables.table('TC00.MS', readonly=False)

print(t.getcoldesc('TEST'))
print(t.getdminfo('TEST'))

