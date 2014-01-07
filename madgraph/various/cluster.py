################################################################################
# Copyright (c) 2009 The MadGraph Development team and Contributors             
#
# This file is a part of the MadGraph 5 project, an application which           
# automatically generates Feynman diagrams and matrix elements for arbitrary    
# high-energy processes in the Standard Model and beyond.                       
#
# It is subject to the MadGraph license which should accompany this             
# distribution.                                                                 
#                                                                               
# For more information, please visit: http://madgraph.phys.ucl.ac.be            
#                                                                               
################################################################################
import subprocess
import logging
import os
import time
import re
import glob

logger = logging.getLogger('madgraph.cluster') 

try:
    from madgraph import MadGraph5Error
    import madgraph.various.misc as misc
except Exception, error:
    if __debug__:
        print  str(error)
    from internal import MadGraph5Error
    import internal.misc as misc

pjoin = os.path.join
   
class ClusterManagmentError(MadGraph5Error):
    pass

class NotImplemented(MadGraph5Error):
    pass


multiple_try = misc.multiple_try

def check_interupt(error=KeyboardInterrupt):

    def deco_interupt(f):
        def deco_f_interupt(self, *args, **opt):
            try:
                return f(self, *args, **opt)
            except error:
                self.remove(*args, **opt)
                raise error
        return deco_f_interupt
    return deco_interupt

class Cluster(object):
    """Basic Class for all cluster type submission"""
    name = 'mother class'

    def __init__(self,*args, **opts):
        """Init the cluster"""

        self.submitted = 0
        self.submitted_ids = []
        self.finish = 0
        
        self.cluster_queue = opts['cluster_queue']
        self.temp_dir = opts['cluster_temp_path']
        self.options = {'cluster_status_update': (600, 30)}
        for key,value in opts.items():
            self.options[key] = value
        

    def submit(self, prog, argument=[], cwd=None, stdout=None, stderr=None, log=None):
        """How to make one submission. Return status id on the cluster."""
        raise NotImplemented, 'No implementation of how to submit a job to cluster \'%s\'' % self.name

    def submit2(self, prog, argument=[], cwd=None, stdout=None, stderr=None, 
                log=None, input_files=[], output_files=[]):
        """How to make one submission. Return status id on the cluster.
        NO SHARE DISK"""

        if cwd is None:
            cwd = os.getcwd()
        if not os.path.exists(prog):
            prog = os.path.join(cwd, prog)
        
        if not hasattr(self, 'temp_dir') or not self.temp_dir:
            return self.submit(prog, argument, cwd, stdout, stderr, log)

        if not input_files and not output_files:
            misc.sprint('not using submit2: no input/output')
            return self.submit(prog, argument, cwd, stdout, stderr, log)

        if cwd is None:
            cwd = os.getcwd()
        if not os.path.exists(prog):
            prog = os.path.join(cwd, prog)
        temp_file_name = "sub." + os.path.basename(prog) + '.'.join(argument)

        text = """#!/bin/bash
        MYTMP=%(tmpdir)s/run$%(job_id)s
        MYPWD=%(cwd)s
        mkdir -p $MYTMP
        cd $MYPWD
        input_files=( %(input_files)s )
        for i in ${input_files[@]}
        do
            cp -R -L $i $MYTMP
        done
        cd $MYTMP
        bash ./%(script)s %(arguments)s
        output_files=( %(output_files)s )
        for i in ${output_files[@]}
        do
            cp -r $MYTMP/$i $MYPWD
        done
        rm -rf $MYTMP
        """
        
        dico = {'tmpdir' : self.temp_dir, 'script': os.path.basename(prog),
                'cwd': cwd, 'job_id': self.job_id,
                'input_files': ' '.join(input_files + [prog]),
                'output_files': ' '.join(output_files),
                'arguments': ' '.join(argument)}
        
        # writing a new script for the submission
        new_prog = pjoin(cwd, temp_file_name)
        open(new_prog, 'w').write(text % dico)
        misc.Popen(['chmod','+x',new_prog],cwd=cwd)
        
        return self.submit(new_prog, argument, cwd, stdout, stderr, log)
        

    def control(self, me_dir=None):
        """Check the status of job associated to directory me_dir. return (idle, run, finish, fail)"""
        if not self.submitted_ids:
            raise NotImplemented, 'No implementation of how to control the job status to cluster \'%s\'' % self.name
        idle, run, fail = 0, 0, 0
        for pid in self.submitted_ids[:]:
            status = self.control_one_job(id)
            if status == 'I':
                idle += 1
            elif status == 'R':
                run += 1
            elif status == 'F':
                self.finish +=1
                self.submitted_ids.remove(pid)
            else:
                fail += 1

        return idle, run, self.finish, fail

    def control_one_job(self, pid):
        """ control the status of a single job with it's cluster id """
        raise NotImplemented, 'No implementation of how to control the job status to cluster \'%s\'' % self.name

    @check_interupt()
    def wait(self, me_dir, fct):
        """Wait that all job are finish"""
        
        nb_iter = 0
        change_at = 5 # number of iteration from which we wait longer between update.
        while 1: 
            nb_iter += 1
            idle, run, finish, fail = self.control(me_dir)
            if fail:
                raise ClusterManagmentError('Some Jobs are in a Hold/... state. Please try to investigate or contact the IT team')
            if idle + run == 0:
                time.sleep(20) #security to ensure that the file are really written on the disk
                logger.info('All jobs finished')
                break
            fct(idle, run, finish)
            if idle < run or nb_iter < change_at:
                time.sleep(self.options['cluster_status_update'][1])
            elif nb_iter == change_at:
                logger.info('''Start to wait %ss between checking status.
Note that you can change this time in the configuration file.
Press ctrl-C to force the update.''' % self.options['cluster_status_update'][0])
                try:
                    time.sleep(self.options['cluster_status_update'][0])
                except KeyboardInterrupt:
                    logger.info('start to update the status')
                    nb_iter = min(0, change_at -2)
            else:
                try:
                    time.sleep(self.options['cluster_status_update'][0])
                except KeyboardInterrupt:
                    logger.info('start to update the status')
                    nb_iter = min(0, change_at -2)
                    
                    
        self.submitted = 0
        self.submitted_ids = []

    @check_interupt()
    def launch_and_wait(self, prog, argument=[], cwd=None, stdout=None, 
                                                         stderr=None, log=None):
        """launch one job on the cluster and wait for it"""
        
        special_output = False # tag for concatenate the error with the output.
        if stderr == -2 and stdout: 
            #We are suppose to send the output to stdout
            special_output = True
            stderr = stdout + '.err'
        pid = self.submit(prog, argument, cwd, stdout, stderr, log)
        while 1:        
            status = self.control_one_job(pid)
            if not status in ['R','I']:
                time.sleep(30) #security to ensure that the file are really written on the disk
                break
            time.sleep(self.options['cluster_status_update'][1])
        
        if special_output:
            # combine the stdout and the stderr
            #wait up to 50 s to see if those files exists
            for i in range(5):
                if os.path.exists(stdout):
                    if not os.path.exists(stderr):
                        time.sleep(5)
                    if os.path.exists(stderr):
                        err_text = open(stderr).read()
                        if not err_text:
                            return
                        logger.warning(err_text)                        
                        text = open(stdout).read()
                        open(stdout,'w').write(text + err_text)
                    else:
                        return
                time.sleep(10)
                        
    def remove(self, *args):
        """ """
        logger.warning("""This cluster didn't support job removal, 
    the jobs are still running on the cluster.""")

class MultiCore(Cluster):
    """ class for dealing with the submission in multiple node"""
    
    job_id = '$'
    
    def __init__(self, *args, **opt):
        """Init the cluster"""
        import thread
        super(MultiCore, self).__init__(self, *args, **opt)
        
        
        self.submitted = 0
        self.finish = 0
        self.nb_core = opt['nb_core']
        self.update_fct = None

        # initialize the thread controler
        self.need_waiting = False
        self.nb_used = 0
        self.lock = thread.allocate_lock()
        self.done = 0 
        self.waiting_submission = []
        self.pids = []
        self.fail_msg = None
        
    def launch_and_wait(self, prog, argument=[], cwd=None, stdout=None, 
                                                         stderr=None, log=None):
        """launch one job and wait for it"""    
        if isinstance(stdout, str):
            stdout = open(stdout, 'w')
        if isinstance(stderr, str):
            stdout = open(stderr, 'w')        
        return misc.call([prog] + argument, stdout=stdout, stderr=stderr, cwd=cwd) 
    
    
    def submit(self, prog, argument=[], cwd=None, stdout=None, stderr=None, 
               log=None):
        """submit a job on multicore machine"""
        
        self.submitted +=1
        if cwd is None:
            cwd = os.getcwd()
        if not os.path.exists(prog) and not misc.which(prog):
            prog = os.path.join(cwd, prog)
        
        import thread
        if self.waiting_submission or self.nb_used == self.nb_core:
            self.waiting_submission.append((prog, argument,cwd, stdout))
            # check that none submission is already finished
            while self.nb_used <  self.nb_core and self.waiting_submission:
                arg = self.waiting_submission.pop(0)
                self.nb_used += 1 # udpate the number of running thread
                thread.start_new_thread(self.launch, arg)              
        elif self.nb_used <  self.nb_core -1:
            self.nb_used += 1 # upate the number of running thread
            thread.start_new_thread(self.launch, (prog, argument, cwd, stdout))
        elif self.nb_used ==  self.nb_core -1:
            self.nb_used += 1 # upate the number of running thread            
            thread.start_new_thread(self.launch, (prog, argument, cwd, stdout))
        
        
    def launch(self, exe, argument, cwd, stdout):
        """ way to launch for multicore."""
        import thread
        def end(self, pid):
            self.nb_used -= 1
            self.done += 1
            try:
                self.pids.remove(pid)
            except:
                pass
            
        fail_msg = None
        try:  
            if os.path.exists(exe) and not exe.startswith('/'):
                exe = './' + exe
            proc = misc.Popen([exe] + argument, cwd=cwd, stdout=stdout, 
                                                           stderr=subprocess.STDOUT)
            pid = proc.pid
            self.pids.append(pid)
            proc.wait()
            if proc.returncode not in [0, 143, -15]:
                fail_msg = 'program %s launch ends with non zero status: %s. Stop all computation' % \
                        (' '.join([exe]+argument), proc.returncode)
                #self.fail_msg = fail_msg
                logger.warning(fail_msg)
                try:
                    log = open(glob.glob(pjoin(cwd,'*','log.txt'))[0]).read()
                    logger.warning('Last 15 lines of lofgile %s:\n%s\n' % \
                            (pjoin(cwd,'*','log.txt'), '\n'.join(log.split('\n')[-15:-1]) + '\n'))
                except IOError, AttributeError:
                    logger.warning('Please look for possible logfiles in %s' % cwd)
                    pass
                self.remove(fail_msg)
            
            # release the lock for allowing to launch the next job
            security = 0       
            # check that the status is locked to avoid coincidence unlock
            while 1:
                while not self.lock.locked():
                    if not self.need_waiting:
                        # Main is not yet locked
                        end(self, pid)
                        return
                    elif security > 60:
                        end(self, pid)
                        return 
                    security += 1
                    time.sleep(1)
                try:
                    self.lock.release()
                except thread.error:
                    continue
                break
            end(self, pid)


        except Exception, error:
            #logger.critical('one core fails with %s' % error)
            self.remove()
            raise

            
          

    def wait(self, me_dir, update_status):
        """Wait that all thread finish
        self.nb_used and self.done are update via each jobs (thread and local)
        self.submitted is the nb of times that submitted has been call (local)
        remaining is the nb of job that we still have to wait. (local)
        self.pids is the list of the BASH pid of the submitted jobs. (thread)
        
        WARNING: In principle all those value are coherent but since some are
        modified in various thread, those data can be corrupted. (not the local 
        one). Nb_used in particular shouldn't be trusted too much.
        This code check in different ways that all jobs have finished.

        In principle, the statement related to  '#security #X' are not used.
        In practise they are times to times.
        """
        
        import thread

        remaining = self.submitted - self.done

        while self.nb_used < self.nb_core:
            if self.waiting_submission:
                arg = self.waiting_submission.pop(0)
                thread.start_new_thread(self.launch, arg)
                self.nb_used += 1 # update the number of running thread
            else:
                break
                    
        try:            
            self.need_waiting = True
            self.lock.acquire()
            no_in_queue = 0
            secure_mode = False # forbid final acauire if in securemode
            while self.waiting_submission or self.nb_used:
                if self.fail_msg:
                    msg,  self.fail_msg = self.fail_msg, None
                    self.remove()
                    raise Exception, msg
                if update_status:
                    update_status(len(self.waiting_submission), self.nb_used, self.done)
                # security#1 that all job expected to be launched since 
                # we enter in this function are indeed launched.
                if len(self.waiting_submission) == 0 == remaining :
                    self.done = self.submitted
                    break
                
                # security #2: nb_used >0 but nothing remains as BASH PID
                if len(self.waiting_submission) == 0 and len(self.pids) == 0:
                    if self.submitted == self.done:
                        break
                    logger.debug('Found too many jobs. Recovering')
                    no_in_queue += 1
                    time.sleep(min(180, 5 * no_in_queue))
                    if no_in_queue > 3:
                        logger.debug('Still too many jobs. Continue')
                        break
                    continue
                
                # security #3: if nb_used not reliable pass in secure mode
                if not secure_mode and len(self.waiting_submission) != 0:
                    if self.nb_used != self.nb_core:
                        if self.nb_used != len(self.pids):
                            secure_mode = True
                # security #4: nb_used not reliable use secure mode to finish the run
                if secure_mode and not self.waiting_submission:
                    self.need_waiting = False
                    if self.lock.locked():
                        self.lock.release()
                    break
                
                # Wait for core to finish               
                self.lock.acquire()
                remaining -=1    # update remaining job
                #submit next one
                if self.waiting_submission:
                    arg = self.waiting_submission.pop(0)
                    thread.start_new_thread(self.launch, arg)
                    self.nb_used += 1 # update the number of running thread

            if self.fail_msg:
                msg,  self.fail_msg = self.fail_msg, None
                self.remove()
                raise Exception, msg            
            # security #5: checked that self.nb_used is not lower than expected
            #This is the most current problem.
            no_in_queue = 0
            while self.submitted > self.done:
                if self.fail_msg:
                    msg,  self.fail_msg = self.fail_msg, None
                    self.remove()
                    raise Exception, msg
                if no_in_queue == 0:
                    logger.debug('Some jobs have been lost. Try to recover')
                #something bad happens
                if not len(self.pids):
                    # The job is not running 
                    logger.critical('Some jobs have been lost in the multicore treatment.')
                    logger.critical('The results might be incomplete. (Trying to continue anyway)')
                    break
                elif update_status:
                    update_status(len(self.waiting_submission), len(self.pids) ,
                                                                      self.done)
                # waiting that those jobs ends.
                if not secure_mode:
                    self.lock.acquire()
                else:
                    no_in_queue += 1
                    try:
                        time.sleep(min(180,5*no_in_queue))
                        if no_in_queue > 5 * 3600.0 / 162:
                            break
                    except KeyboardInterrupt:
                        logger.warning('CTRL-C assumes that all jobs are done. Continue the code')
                        self.pids = [] # avoid security 6
                        break
                    
            # security #6. check that queue is empty. don't
            no_in_queue = 0
            while len(self.pids):
                if self.fail_msg:
                    msg,  self.fail_msg = self.fail_msg, None
                    self.remove()
                    raise Exception, msg
                self.need_waiting = False
                if self.lock.locked():
                        self.lock.release()
                secure_mode = True
                if no_in_queue == 0 : 
                    logger.warning('Some jobs have been lost. Try to recover.')
                    logger.warning('Hitting ctrl-c will consider that all jobs are done and continue the code.')
                try:
                    #something very bad happens
                    if update_status:
                        update_status(len(self.waiting_submission), len(self.pids) ,
                                                                      self.done)
                    time.sleep(min(5*no_in_queue, 180))
                    no_in_queue += 1
                    if no_in_queue > 5 * 3600.0 / 162:
                            break
                except KeyboardInterrupt:
                    break
                
            # print a last time the status (forcing 0 for the running)  
            if update_status:
                self.next_update = 0
                update_status(len(self.waiting_submission), 0, self.done)             
            
            # reset variable for next submission
            self.need_waiting = False
            security = 0 
            while not self.lock.locked() and security < 10:
                # check that the status is locked to avoid coincidence unlock
                if secure_mode:
                    security = 10
                security +=1
                time.sleep(1)
            if security < 10:
                self.lock.release()
            self.done = 0
            self.nb_used = 0
            self.submitted = 0
            self.pids = []
            
        except KeyboardInterrupt:
            self.remove()
            raise
        if self.fail_msg:
            msg,  self.fail_msg = self.fail_msg, None
            self.remove()
            raise Exception, msg 
            
    def remove(self, error=None):
        """Ensure that all thread are killed"""
        logger.info('remove job currently running')
        self.waiting_submission = []
        if error:
            self.fail_msg = error
        for pid in list(self.pids):
            out = os.system('CPIDS=$(pgrep -P %(pid)s); kill -15 $CPIDS > /dev/null 2>&1' \
                            % {'pid':pid} )
            out = os.system('kill -15 %(pid)s > /dev/null 2>&1' % {'pid':pid} )            
            if out == 0:
                try:
                    self.pids.remove(pid)
                except:
                    pass
            #out = os.system('kill -9 %s &> /dev/null' % pid)

        time.sleep(1) # waiting if some were submitting at the time of ctrl-c
        for pid in list(self.pids):
            out = os.system('CPIDS=$(pgrep -P %s); kill -15 $CPIDS > /dev/null 2>&1' % pid )
            out = os.system('kill -15 %(pid)s > /dev/null 2>&1' % {'pid':pid} ) 
            if out == 0:
                try:
                    self.pids.remove(pid)
                except:
                    pass
                    
class CondorCluster(Cluster):
    """Basic class for dealing with cluster submission"""
    
    name = 'condor'
    job_id = 'CONDOR_ID'

    @multiple_try()
    def submit(self, prog, argument=[], cwd=None, stdout=None, stderr=None, log=None):
        """Submit a job prog to a Condor cluster"""
        
        text = """Executable = %(prog)s
                  output = %(stdout)s
                  error = %(stderr)s
                  log = %(log)s
                  %(argument)s
                  environment = CONDOR_ID=$(Cluster).$(Process)
                  Universe = vanilla
                  notification = Error
                  Initialdir = %(cwd)s
                  %(requirement)s
                  getenv=True
                  queue 1
               """
        
        if self.cluster_queue not in ['None', None]:
            requirement = 'Requirements = %s=?=True' % self.cluster_queue
        else:
            requirement = ''

        if cwd is None:
            cwd = os.getcwd()
        if stdout is None:
            stdout = '/dev/null'
        if stderr is None:
            stderr = '/dev/null'
        if log is None:
            log = '/dev/null'
        if not os.path.exists(prog):
            prog = os.path.join(cwd, prog)
        if argument:
            argument = 'Arguments = %s' % ' '.join(argument)
        else:
            argument = ''
        

        dico = {'prog': prog, 'cwd': cwd, 'stdout': stdout, 
                'stderr': stderr,'log': log,'argument': argument,
                'requirement': requirement}

        open('submit_condor','w').write(text % dico)
        a = misc.Popen(['condor_submit','submit_condor'], stdout=subprocess.PIPE)
        output = a.stdout.read()
        #Submitting job(s).
        #Logging submit event(s).
        #1 job(s) submitted to cluster 2253622.
        pat = re.compile("submitted to cluster (\d*)",re.MULTILINE)
        try:
            id = pat.search(output).groups()[0]
        except:
            raise ClusterManagmentError, 'fail to submit to the cluster: \n%s' \
                                                                        % output 
        self.submitted += 1
        self.submitted_ids.append(id)
        return id

    @multiple_try()
    def submit2(self, prog, argument=[], cwd=None, stdout=None, stderr=None, 
                log=None, input_files=[], output_files=[]):
        """Submit the job on the cluster NO SHARE DISK
           input/output file should be give relative to cwd
        """
        
        text = """Executable = %(prog)s
                  output = %(stdout)s
                  error = %(stderr)s
                  log = %(log)s
                  %(argument)s
                  should_transfer_files = YES
                  when_to_transfer_output = ON_EXIT
                  transfer_input_files = %(input_files)s
                  %(output_files)s
                  Universe = vanilla
                  notification = Error
                  Initialdir = %(cwd)s
                  %(requirement)s
                  getenv=True
                  queue 1
               """
        
        if self.cluster_queue not in ['None', None]:
            requirement = 'Requirements = %s=?=True' % self.cluster_queue
        else:
            requirement = ''

        if cwd is None:
            cwd = os.getcwd()
        if stdout is None:
            stdout = '/dev/null'
        if stderr is None:
            stderr = '/dev/null'
        if log is None:
            log = '/dev/null'
        if not os.path.exists(prog):
            prog = os.path.join(cwd, prog)
        if argument:
            argument = 'Arguments = %s' % ' '.join(argument)
        else:
            argument = ''
        # input/output file treatment
        if input_files:
            input_files = ','.join(input_files)
        else: 
            input_files = ''
        if output_files:
            output_files = 'transfer_output_files = %s' % ','.join(output_files)
        else:
            output_files = ''
        
        

        dico = {'prog': prog, 'cwd': cwd, 'stdout': stdout, 
                'stderr': stderr,'log': log,'argument': argument,
                'requirement': requirement, 'input_files':input_files, 
                'output_files':output_files}

        open('submit_condor','w').write(text % dico)
        a = subprocess.Popen(['condor_submit','submit_condor'], stdout=subprocess.PIPE)
        output = a.stdout.read()
        #Submitting job(s).
        #Logging submit event(s).
        #1 job(s) submitted to cluster 2253622.
        pat = re.compile("submitted to cluster (\d*)",re.MULTILINE)
        try:
            id = pat.search(output).groups()[0]
        except:
            raise ClusterManagmentError, 'fail to submit to the cluster: \n%s' \
                                                                        % output 
        self.submitted += 1
        self.submitted_ids.append(id)
        return id




    
    @multiple_try(nb_try=10, sleep=10)
    def control_one_job(self, id):
        """ control the status of a single job with it's cluster id """
        cmd = 'condor_q '+str(id)+" -format \'%-2s \\n\' \'ifThenElse(JobStatus==0,\"U\",ifThenElse(JobStatus==1,\"I\",ifThenElse(JobStatus==2,\"R\",ifThenElse(JobStatus==3,\"X\",ifThenElse(JobStatus==4,\"C\",ifThenElse(JobStatus==5,\"H\",ifThenElse(JobStatus==6,\"E\",string(JobStatus))))))))\'"
        status = misc.Popen([cmd], shell=True, stdout=subprocess.PIPE, 
                                                         stderr=subprocess.PIPE)
        
        error = status.stderr.read()
        if status.returncode or error:
            raise ClusterManagmentError, 'condor_q returns error: %s' % error

        return status.stdout.readline().strip()
    
    @check_interupt()
    @multiple_try(nb_try=10, sleep=10)
    def control(self, me_dir):
        """ control the status of a single job with it's cluster id """
        
        if not self.submitted_ids:
            return 0, 0, 0, 0
        
        packet = 15000
        for i in range(1+(len(self.submitted_ids)-1)//packet):
            start = i * packet
            stop = (i+1) * packet
            cmd = "condor_q " + ' '.join(self.submitted_ids[start:stop]) + " -format \'%-2s \\n\' \'ifThenElse(JobStatus==0,\"U\",ifThenElse(JobStatus==1,\"I\",ifThenElse(JobStatus==2,\"R\",ifThenElse(JobStatus==3,\"X\",ifThenElse(JobStatus==4,\"C\",ifThenElse(JobStatus==5,\"H\",ifThenElse(JobStatus==6,\"E\",string(JobStatus))))))))\'"
            
                
            
            status = misc.Popen([cmd], shell=True, stdout=subprocess.PIPE, 
                                                             stderr=subprocess.PIPE)
            error = status.stderr.read()
            if status.returncode or error:
                raise ClusterManagmentError, 'condor_q returns error: %s' % error
                
                
            idle, run, fail = 0, 0, 0
            for line in status.stdout:
                status = line.strip()
                if status in ['I','U']:
                    idle += 1
                elif status == 'R':
                    run += 1
                elif status != 'C':
                    fail += 1

        return idle, run, self.submitted - (idle+run+fail), fail
    
    @multiple_try()
    def remove(self, *args):
        """Clean the jobson the cluster"""
        
        if not self.submitted_ids:
            return
        cmd = "condor_rm %s" % ' '.join(self.submitted_ids)
        
        status = misc.Popen([cmd], shell=True, stdout=open(os.devnull,'w'))
        
class PBSCluster(Cluster):
    """Basic class for dealing with cluster submission"""
    
    name = 'pbs'
    job_id = 'PBS_JOBID'
    idle_tag = ['Q']
    running_tag = ['T','E','R']
    complete_tag = ['C']

    @multiple_try()
    def submit(self, prog, argument=[], cwd=None, stdout=None, stderr=None, log=None):
        """Submit a job prog to a PBS cluster"""
        
        me_dir = os.path.realpath(os.path.join(cwd,prog)).rsplit('/SubProcesses',1)[0]
        me_dir = misc.digest(me_dir)[-14:]
        if not me_dir[0].isalpha():
            me_dir = 'a' + me_dir[1:]
        
        text = ""
        if cwd is None:
            cwd = os.getcwd()
        else: 
            text = " cd %s;" % cwd
        if stdout is None:
            stdout = '/dev/null'
        if stderr is None:
            stderr = '/dev/null'
        elif stderr == -2: # -2 is subprocess.STDOUT
            stderr = stdout
        if log is None:
            log = '/dev/null'
        
        if not os.path.isabs(prog):
            text += "./%s" % prog
        else:
            text+= prog
        
        if argument:
            text += ' ' + ' '.join(argument)

        command = ['qsub','-o', stdout,
                   '-N', me_dir, 
                   '-e', stderr,
                   '-V']

        if self.cluster_queue and self.cluster_queue != 'None':
            command.extend(['-q', self.cluster_queue])

        a = misc.Popen(command, stdout=subprocess.PIPE, 
                                      stderr=subprocess.STDOUT,
                                      stdin=subprocess.PIPE, cwd=cwd)
            
        output = a.communicate(text)[0]
        id = output.split('.')[0]
        if not id.isdigit():
            raise ClusterManagmentError, 'fail to submit to the cluster: \n%s' \
                                                                        % output 
        self.submitted += 1
        self.submitted_ids.append(id)
        return id

    @multiple_try()
    def control_one_job(self, id):
        """ control the status of a single job with it's cluster id """
        cmd = 'qstat '+str(id)
        status = misc.Popen([cmd], shell=True, stdout=subprocess.PIPE,
                                  stderr=open(os.devnull,'w'))
        
        for line in status.stdout:
            line = line.strip()
            if 'Unknown' in line:
                return 'F'
            elif line.startswith(str(id)):
                status = line.split()[4]
        if status in self.idle_tag:
            return 'I' 
        elif status in self.running_tag:                
            return 'R' 
        return 'F'
        

    @multiple_try()    
    def control(self, me_dir):
        """ control the status of a single job with it's cluster id """
        cmd = "qstat"
        status = misc.Popen([cmd], stdout=subprocess.PIPE)

        if me_dir.endswith('/'):
            me_dir = me_dir[:-1]    
        me_dir = misc.digest(me_dir)[-14:]
        if not me_dir[0].isalpha():
            me_dir = 'a' + me_dir[1:]

        idle, run, fail = 0, 0, 0
        for line in status.stdout:
            if me_dir in line:
                status = line.split()[4]
                if status in self.idle_tag:
                    idle += 1
                elif status in self.running_tag:
                    run += 1
                elif status in self.complete_tag:
                    continue
                else:
                    fail += 1

        return idle, run, self.submitted - (idle+run+fail), fail

    @multiple_try()
    def remove(self, *args):
        """Clean the jobs on the cluster"""
        
        if not self.submitted_ids:
            return
        cmd = "qdel %s" % ' '.join(self.submitted_ids)
        status = misc.Popen([cmd], shell=True, stdout=open(os.devnull,'w'))


class SGECluster(Cluster):
    """Basic class for dealing with cluster submission"""
    # Class written by Arian Abrahantes.

    name = 'sge'
    job_id = 'JOB_ID'
    idle_tag = ['qw', 'hqw','hRqw','w']
    running_tag = ['r','t','Rr','Rt']

    def def_get_path(self,location):
        """replace string for path issues"""
        location = os.path.realpath(location)
        homePath = os.getenv("HOME")
        if homePath:
            location = location.replace(homePath,'$HOME')
        return location

    @multiple_try()
    def submit(self, prog, argument=[], cwd=None, stdout=None, stderr=None, log=None):
        """Submit a job prog to an SGE cluster"""

        me_dir = os.path.realpath(os.path.join(cwd,prog)).rsplit('/SubProcesses',1)[0]
        me_dir = misc.digest(me_dir)[-10:]
        if not me_dir[0].isalpha():
            me_dir = 'a' + me_dir[1:]

        if cwd is None:
            #cwd = os.getcwd()
            cwd = self.def_get_path(os.getcwd())
        cwd1 = self.def_get_path(cwd)
        text = " cd %s;" % cwd1
        if stdout is None:
            stdout = '/dev/null'
        else:
            stdout = self.def_get_path(stdout)
        if stderr is None:
            stderr = '/dev/null'
        elif stderr == -2: # -2 is subprocess.STDOUT
            stderr = stdout
        else:
            stderr = self.def_get_path(stderr)
            
        if log is None:
            log = '/dev/null'
        else:
            log = self.def_get_path(log)

        text += prog
        if argument:
            text += ' ' + ' '.join(argument)

        #if anything slips through argument
        #print "!=== inteded change ",text.replace('/srv/nfs','')
        #text = text.replace('/srv/nfs','')
        homePath = os.getenv("HOME")
        if homePath:
            text = text.replace(homePath,'$HOME')

        logger.debug("!=== input  %s" % text)
        logger.debug("!=== output %s" %  stdout)
        logger.debug("!=== error  %s" % stderr)
        logger.debug("!=== logs   %s" % log)

        command = ['qsub','-o', stdout,
                   '-N', me_dir, 
                   '-e', stderr,
                   '-V']

        if self.cluster_queue and self.cluster_queue != 'None':
            command.extend(['-q', self.cluster_queue])

        a = misc.Popen(command, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             stdin=subprocess.PIPE, cwd=cwd)

        output = a.communicate(text)[0]
        id = output.split(' ')[2]
        if not id.isdigit():
            raise ClusterManagmentError, 'fail to submit to the cluster: \n%s' \
                                                                        % output 
        self.submitted += 1
        self.submitted_ids.append(id)
        logger.debug(output)

        return id

    @multiple_try()
    def control_one_job(self, id):
        """ control the status of a single job with it's cluster id """
        #cmd = 'qstat '+str(id)
        cmd = 'qstat '
        status = misc.Popen([cmd], shell=True, stdout=subprocess.PIPE)
        for line in status.stdout:
            #print "!==",line
            #line = line.strip()
            #if 'Unknown' in line:
            #    return 'F'
            #elif line.startswith(str(id)):
            #    status = line.split()[4]
            if str(id) in line:
                status = line.split()[4]
                #print "!=status", status
        if status in self.idle_tag:
            return 'I' 
        elif status in self.running_tag:                
            return 'R' 
        return 'F'

    @multiple_try()
    def control(self, me_dir):
        """ control the status of a single job with it's cluster id """
        cmd = "qstat "
        status = misc.Popen([cmd], shell=True, stdout=subprocess.PIPE)

        if me_dir.endswith('/'):
           me_dir = me_dir[:-1]    
        me_dir = misc.digest(me_dir)[-10:]
        if not me_dir[0].isalpha():
            me_dir = 'a' + me_dir[1:]

        idle, run, fail = 0, 0, 0
        for line in status.stdout:
            if me_dir in line:
                status = line.split()[4]
                if status in self.idle_tag:
                    idle += 1
                elif status in self.running_tag:
                    run += 1
                else:
                    logger.debug(line)
                    fail += 1

        return idle, run, self.submitted - (idle+run+fail), fail

    
    
    @multiple_try()
    def remove(self, *args):
        """Clean the jobs on the cluster"""
        
        if not self.submitted_ids:
            return
        cmd = "qdel %s" % ' '.join(self.submitted_ids)
        status = misc.Popen([cmd], shell=True, stdout=open(os.devnull,'w'))


class LSFCluster(Cluster):
    """Basic class for dealing with cluster submission"""
    
    name = 'lsf'
    job_id = 'LSB_JOBID'

    @multiple_try()
    def submit(self, prog, argument=[], cwd=None, stdout=None, stderr=None, log=None):
        """Submit the job prog to an LSF cluster"""
        
        me_dir = os.path.realpath(os.path.join(cwd,prog)).rsplit('/SubProcesses',1)[0]
        me_dir = misc.digest(me_dir)[-14:]
        if not me_dir[0].isalpha():
            me_dir = 'a' + me_dir[1:]
        
        text = ""
        command = ['bsub', '-J', me_dir]
        if cwd is None:
            cwd = os.getcwd()
        else: 
            text = " cd %s;" % cwd
        if stdout and isinstance(stdout, str):
            command.extend(['-o', stdout])
        if stderr and isinstance(stdout, str):
            command.extend(['-e', stderr])
        elif stderr == -2: # -2 is subprocess.STDOUT
            pass
        if log is None:
            log = '/dev/null'
        
        text += prog
        if argument:
            text += ' ' + ' '.join(argument)
        
        if self.cluster_queue and self.cluster_queue != 'None':
            command.extend(['-q', self.cluster_queue])

        a = misc.Popen(command, stdout=subprocess.PIPE, 
                                      stderr=subprocess.STDOUT,
                                      stdin=subprocess.PIPE, cwd=cwd)
            
        output = a.communicate(text)[0]
        #Job <nnnn> is submitted to default queue <normal>.
        try:
            id = output.split('>',1)[0].split('<')[1]
        except:
            raise ClusterManagmentError, 'fail to submit to the cluster: \n%s' \
                                                                        % output 
        if not id.isdigit():
            raise ClusterManagmentError, 'fail to submit to the cluster: \n%s' \
                                                                        % output 
        self.submitted += 1
        self.submitted_ids.append(id)
        return id        
        
        
    @multiple_try()
    def control_one_job(self, id):
        """ control the status of a single job with it's cluster id """
        
        cmd = 'bjobs '+str(id)
        status = misc.Popen([cmd], shell=True, stdout=subprocess.PIPE)
        
        for line in status.stdout:
            line = line.strip().upper()
            if 'JOBID' in line:
                continue
            elif str(id) not in line:
                continue
            status = line.split()[2]
            if status == 'RUN':
                return 'R'
            elif status == 'PEND':
                return 'I'
            elif status == 'DONE':
                return 'F'
            else:
                return 'H'
            return 'F'

    @multiple_try()   
    def control(self, me_dir):
        """ control the status of a single job with it's cluster id """
        
        if not self.submitted_ids:
            return 0, 0, 0, 0
        
        cmd = "bjobs " + ' '.join(self.submitted_ids) 
        status = misc.Popen([cmd], shell=True, stdout=subprocess.PIPE)

        idle, run, fail = 0, 0, 0
        for line in status.stdout:
            line = line.strip()
            if 'JOBID' in line:
                continue
            splitline = line.split()
            id = splitline[0]
            if id not in self.submitted_ids:
                continue
            status = splitline[2]
            if status == 'RUN':
                run += 1
            elif status == 'PEND':
                idle += 1
            elif status == 'DONE':
                self.submitted_ids.remove(id)
            else:
                fail += 1

        return idle, run, self.submitted - (idle+run+fail), fail

    @multiple_try()
    def remove(self, *args):
        """Clean the jobs on the cluster"""
        
        if not self.submitted_ids:
            return
        cmd = "bkill %s" % ' '.join(self.submitted_ids)
        status = misc.Popen([cmd], shell=True, stdout=open(os.devnull,'w'))

class GECluster(Cluster):
    """Class for dealing with cluster submission on a GE cluster"""
    
    name = 'ge'
    job_id = 'JOB_ID'
    idle_tag = ['qw']
    running_tag = ['r']

    @multiple_try()
    def submit(self, prog, argument=[], cwd=None, stdout=None, stderr=None, log=None):
        """Submit a job prog to a GE cluster"""
        
        text = ""
        if cwd is None:
            cwd = os.getcwd()
        else: 
            text = " cd %s; bash " % cwd
        if stdout is None:
            stdout = os.path.join(cwd, "log.%s" % prog.split('/')[-1])
        if stderr is None:
            stderr = os.path.join(cwd, "err.%s" % prog.split('/')[-1])
        elif stderr == -2: # -2 is subprocess.STDOUT
            stderr = stdout
        if log is None:
            log = '/dev/null'

        text += prog
        if argument:
            text += ' ' + ' '.join(argument)
        text += '\n'
        tmp_submit = os.path.join(cwd, 'tmp_submit')
        open(tmp_submit,'w').write(text)

        a = misc.Popen(['qsub','-o', stdout,
                                     '-e', stderr,
                                     tmp_submit],
                                     stdout=subprocess.PIPE, 
                                     stderr=subprocess.STDOUT,
                                     stdin=subprocess.PIPE, cwd=cwd)

        output = a.communicate()[0]
        #Your job 874511 ("test.sh") has been submitted
        pat = re.compile("Your job (\d*) \(",re.MULTILINE)
        try:
            id = pat.search(output).groups()[0]
        except:
            raise ClusterManagmentError, 'fail to submit to the cluster: \n%s' \
                                                                        % output 
        self.submitted += 1
        self.submitted_ids.append(id)
        return id

    @multiple_try()
    def control_one_job(self, id):
        """ control the status of a single job with it's cluster id """
        cmd = 'qstat | grep '+str(id)
        status = misc.Popen([cmd], shell=True, stdout=subprocess.PIPE)
        if not status:
            return 'F'
        #874516 0.00000 test.sh    alwall       qw    03/04/2012 22:30:35                                    1
        pat = re.compile("^(\d+)\s+[\d\.]+\s+[\w\d\.]+\s+[\w\d\.]+\s+(\w+)\s")
        stat = ''
        for line in status.stdout.read().split('\n'):
            if not line:
                continue
            line = line.strip()
            try:
                groups = pat.search(line).groups()
            except:
                raise ClusterManagmentError, 'bad syntax for stat: \n\"%s\"' % line
            if groups[0] != id: continue
            stat = groups[1]
        if not stat:
            return 'F'
        if stat in self.idle_tag:
            return 'I' 
        if stat in self.running_tag:                
            return 'R' 
        
    @multiple_try()
    def control(self, me_dir=None):
        """Check the status of job associated to directory me_dir. return (idle, run, finish, fail)"""
        if not self.submitted_ids:
            return 0, 0, 0, 0
        idle, run, fail = 0, 0, 0
        ongoing = []
        for statusflag in ['p', 'r', 'sh']:
            cmd = 'qstat -s %s' % statusflag
            status = misc.Popen([cmd], shell=True, stdout=subprocess.PIPE)
            #874516 0.00000 test.sh    alwall       qw    03/04/2012 22:30:35                                    1
            pat = re.compile("^(\d+)")
            for line in status.stdout.read().split('\n'):
                line = line.strip()
                try:
                    id = pat.search(line).groups()[0]
                except Exception:
                    pass
                else:
                    if id not in self.submitted_ids:
                        continue
                    ongoing.append(id)
                    if statusflag == 'p':
                        idle += 1
                    if statusflag == 'r':
                        run += 1
                    if statusflag == 'sh':
                        fail += 1

        self.submitted_ids = ongoing

        return idle, run, self.submitted - idle - run - fail, fail

    @multiple_try()
    def remove(self, *args):
        """Clean the jobs on the cluster"""
        
        if not self.submitted_ids:
            return
        cmd = "qdel %s" % ' '.join(self.submitted_ids)
        status = misc.Popen([cmd], shell=True, stdout=open(os.devnull,'w'))

def asyncrone_launch(exe, cwd=None, stdout=None, argument = [], **opt):
    """start a computation and not wait for it to finish.
       this fonction returns a lock which is locked as long as the job is 
       running."""

    mc = MultiCore(1)
    mc.submit(exe, argument, cwd, stdout, **opt)
    mc.need_waiting = True
    mc.lock.acquire()
    return mc.lock


class SLURMCluster(Cluster):
    """Basic class for dealing with cluster submission"""

    name = 'slurm'
    job_id = 'SLURM_JOBID'
    idle_tag = ['Q','PD','S','CF']
    running_tag = ['R']
    complete_tag = ['C']

    @multiple_try()
    def submit(self, prog, argument=[], cwd=None, stdout=None, stderr=None, log=None):
        """Submit a job prog to a SLURM cluster"""
        
        me_dir = os.path.realpath(os.path.join(cwd,prog)).rsplit('/SubProcesses',1)[0]
        me_dir = misc.digest(me_dir)[-8:]

        if not me_dir[0].isalpha():
            me_dir = 'a' + me_dir[1:]
        
        if cwd is None:
            cwd = os.getcwd()
        if stdout is None:
            stdout = '/dev/null'
        if stderr is None:
            stderr = '/dev/null'
        elif stderr == -2: # -2 is subprocess.STDOUT
            stderr = stdout
        if log is None:
            log = '/dev/null'
        
        command = ['sbatch','-o', stdout,
                   '-J', me_dir, 
                   '-e', stderr, prog]
                   
        a = misc.Popen(command, stdout=subprocess.PIPE, 
                                      stderr=subprocess.STDOUT,
                                      stdin=subprocess.PIPE, cwd=cwd)

        output = a.communicate()
        output_arr = output[0].split(' ')
        id = output_arr[3].rstrip()

        if not id.isdigit():
            raise ClusterManagmentError, 'fail to submit to the cluster: \n%s' \

        self.submitted += 1
        self.submitted_ids.append(id)
        return id

    @multiple_try()
    def control_one_job(self, id):
        """ control the status of a single job with it's cluster id """
        cmd = 'squeue j'+str(id)
        status = misc.Popen([cmd], shell=True, stdout=subprocess.PIPE,
                                  stderr=open(os.devnull,'w'))
        
        for line in status.stdout:
            line = line.strip()
            if 'Invalid' in line:
                return 'F'
            elif line.startswith(str(id)):
                status = line.split()[4]
        if status in self.idle_tag:
            return 'I' 
        elif status in self.running_tag:                
            return 'R' 
        return 'F'
        
    @multiple_try()    
    def control(self, me_dir):
        """ control the status of a single job with it's cluster id """
        cmd = "squeue"
        status = misc.Popen([cmd], stdout=subprocess.PIPE)

        if me_dir.endswith('/'):
           me_dir = me_dir[:-1]   
        me_dir = misc.digest(me_dir)[-8:]
        if not me_dir[0].isalpha():
                  me_dir = 'a' + me_dir[1:]

        idle, run, fail = 0, 0, 0
        for line in status.stdout:
            if me_dir in line:
                status = line.split()[4]
                if status in self.idle_tag:
                    idle += 1
                elif status in self.running_tag:
                    run += 1
                elif status in self.complete_tag:
                    continue
                else:
                    fail += 1

        return idle, run, self.submitted - (idle+run+fail), fail

    @multiple_try()
    def remove(self, *args):
        """Clean the jobs on the cluster"""
        
        if not self.submitted_ids:
            return
        cmd = "scancel %s" % ' '.join(self.submitted_ids)
        print 'cmd = ',cmd
        status = misc.Popen([cmd], shell=True, stdout=open(os.devnull,'w'))


from_name = {'condor':CondorCluster, 'pbs': PBSCluster, 'sge': SGECluster, 
             'lsf': LSFCluster, 'ge':GECluster, 'slurm': SLURMCluster}


