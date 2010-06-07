# $Id$
# vmd remote control --- client/server scripts to run VMD from python
# Copyright (c) 2007-2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Lesser Public License, version 3 or later.
#
# using asynchat/asyncore to talk to a VMD server process
# (see http://www.python.org/doc/current/lib/module-asyncore.html & 
# http://squirl.nightmare.com/medusa/async_sockets.html)
__docformat__ = "restructuredtext en"

import asynchat, asyncore, socket
import os, sys, readline, string, re

try:
    # part of a python egg
    # see http://peak.telecommunity.com/DevCenter/PythonEggs#accessing-package-resources
    from pkg_resources import resource_filename
    REMOTE_TCL = resource_filename(__name__,'remote_ctl.tcl')
except ImportError:
    REMOTE_TCL = os.path.join(os.path.split(__file__)[0],'remote_ctl.tcl')

services = { 'vmd' : 5555 }
RESPONSE_TERMINATOR = '__END_OF_VMD_RESPONSE__\r\n'  # see remote_ctl.tcl


class server:
    """The VMD server process."""
    def __init__(self,vmdbinary='vmd',
                 server_tcl=REMOTE_TCL,
                 force=False,maxdelay=10,dispdev='text'):
        
        """Start VMD in text mode and launch the remote server.

        :Arguments:
           vmdbinary
                 ``vmd`` (either in PATH or absolute path)
           server_tcl
                 path to the tcl script that starts the listening server in VMD
                 (the default is correct in 99.9% of cases)
           force
                 ``False``
                       don't launch new server if one is already active
                 ``True``
                       always start new vmd process
           maxdelay         
                 maximum time to wait for the server to come up in seconds
           dispdev
                 VMD display device; default is 'text' which runs VMD without
                 graphics. 'win' is the graphical window device

        :Bugs: Starting multiple VMD processes does not work as we always
               use the same port

        """
        devices = {'graphics':'win','win':'win', 'text':'text','batch':'text'}
        self.vmdbinary = vmdbinary
        try:
            self.dispdev = devices[dispdev]
        except KeyError:
            raise ValueError("VMD display device '"+str(dispdev)+"' is not recognized. "+
                             "Choose one of "+str(devices.keys())+".")
        
        # tcl file is stored in the same directory as this module file
        self.server_tcl = server_tcl
        self.vmdbinary = vmdbinary
        self.startcmd = '%s -dispdev %s -e %s' % (self.vmdbinary, self.dispdev, self.server_tcl)
        self.start(force=force)

    def start(self,force=False,maxdelay=10):
        """Start VMD and launch the remote server.

        :Arguments:
           force
               ``False``
                     don't launch new server if one is already active
               ``True``
                     always start new vmd process
           maxdelay
               maximum number of seconds to wait for VMD to start 
        """
        import time
        
        if force or not self.ping():
           args = self.startcmd.split(" ")
           pid = os.spawnvp(os.P_NOWAIT,args[0],args) 
        # now wait until the server is up (check every 2 seconds)
        t = interval = 2        
        time.sleep(interval)
        while not self.ping() and t < maxdelay:
            time.sleep(interval)
            t += interval
        if not self.ping():
            raise RuntimeError('Failed to bring up the VMD server.')
        
    def stop(self):
        """Shutdown VMD."""
        command('quit')

    def ping(self,pid=os.getpid()):
        """Check if a vmd server can be used.

        :Returns: ``True`` for a live VMD server, or ``False``.

        Ignore the message 'error: uncaptured python exception' if the server is down.
        """
        token = 'ALIVE (ping from pid %d)' % pid
        c = command('set __pingtest__ {%s}' % token)  
        x, = c.results()
        return x == token

    def command(self,*args):
        """Send commands to the VMD server::

           c = command('cd','set w [atomselect top {water}]', '$w writepdb water.pdb')
           c.results()

        This is only a thin convenience wrapper for
        :meth:`vmd.control.command` and not strongly tied to the server (as
        anyone can connect).
        """
        return command(*args)

class client(asynchat.async_chat):
    """one command -> response exchange between the client and vmd::

       c = client(host,port=5555)
       c.cmd(tcl, tcl,...)
       asyncore.loop()

    Starting VMD as ``vmd -e remote_ctl.tcl`` opens port 5555 for connection.
    The client only becomes active in the ``asyncore.loop()`` and exits after
    sending the commands and receiving the response. The response is available
    as :meth:`client.response`

    :Parameters:
       host
           currently remote_ctl.tcl only allows 'localhost'
       port
           port to connect to (typically 5555)

    :Methods:
      :meth:`client.cmd`
        commands (with embedded newlines!) scheduled for sending and execution in VMD
      :meth:`client.response`
        response from VMD

    :Bugs: Somehow it doesnt like many commands...
    """
    def __init__(self, host,port=services['vmd']):
        asynchat.async_chat.__init__(self)
        self.port = port
        self.ibuffer = ""
        self.create_socket(socket.AF_INET, socket.SOCK_STREAM)
        self.connect( (host,port) )
        self.set_terminator(RESPONSE_TERMINATOR)

    def handle_connect(self):
        pass
    
    def collect_incoming_data(self,data):
        """buffer incoming data"""
        self.ibuffer += data

    def found_terminator(self):
        self.push('close\n')     # tell the other side we are done
                                 # and have them shut down this socket
                                 # -- THIS IS PART OF THE PROTOCOL

    def cmd(self,*tcl):
        """Submits the commands to be executed in VMD::

           c.cmd(tcl, tcl, ...)    

        Commands (*with embedded newlines!*) scheduled for sending and execution
        in ``VMD``. All strings will be executed sequentially.
        """
        if len(tcl) == 0:
            raise ValueError, 'at least one Tcl command is required'
        s = " ".join(tcl)
        self.push(s)

    def response(self):
        return self.ibuffer
        

class command(client):
    """Send one or more tcl commands to VMD and return response::

       c = command(*tcl)

    Appends a newline to each command if necessary and then
    feeds every single command separately to vmd. Commands that
    include newlines are split on newlines. The responses are stored
    in the tuple c._results (and can be retrieved by the c.results()
    method).

    Technically, this is unelegant cr^&...

    :Methods:
      :meth:`client.results`
         response from VMD
      :meth:`clinet.commands`
         corresponding commands
    """

    def __init__(self,*commands):
        self._resp = []
        commands = [c.strip()+'\n' for c in commands]
        s = " ".join(commands)             # one string to join them all...
        s = re.sub('\n$','',s)             # get rid of last newline        
        self._cmds = string.split(s,'\n')  # ...so that we get only nice command bits
        for c in self._cmds:
            # print "DEBUG: sending command [%s]" % c
            cl = client('localhost')       # what a waste--single command exchange
            cl.cmd(c + "\n")
            asyncore.loop()
            r = cl.response()
            r = re.sub('\r\n','\n',r)
            r = re.sub('\n$','',r)
            # print "DEBUG: response        [%s]" % r
            self._resp.append(r[:])

    def results(self):
        """Return results from vmd."""
        return tuple(self._resp)

    def commands(self):
        """Return submitted tcl commands."""
        return tuple(self._cmds)

class interactive(client):
    """Interactive remote session with vmd::

      interactive(host)
      asyncore.loop()

    When the ``loop()`` is called the interactive session starts and the prompt
    is displayed as ``python->vmd>``. You are now connected to the tcl
    interpreter in vmd. End the session by issuing the command ``close`` or
    ``EXIT``.

    :Parameters:
       host
           currently remote_ctl.tcl only allows 'localhost'    

    Commands interpreted by the remote server and not by tcl in vmd:

    ================  ====================================================
    command           description
    ================  ====================================================
    close             close the current connection (socket)
    exit              exit the server (no more connections possible,
                      but current connections are still open)
    loglevel N        set LOGLEVEL to value N (0<=N<=2) [default: 1]
    EXIT              exit the interactive session
    ================  ====================================================
    """
    
    def __init__(self, host,port=services['vmd']):
        client.__init__(self,host)
        self.stack_count = 0         # number of command responses to expect

    def handle_connect(self):
        # we need one server response
        self.cmd('puts "Interactive connection from client established."')

    def found_terminator(self):
        data = self.ibuffer
        self.ibuffer = ""
        if data.endswith('\r'):
            data = data[:-1]
        print data,
        self.stack_count -= 1

        if self.stack_count == 0:
            # enter command processing if no more output
            self.cmd()
        
    def cmd(self,*commands):
        if len(commands) == 0:
            c = self.getinput()
            if c == 'EXIT\n':
                print "Exit interactive loop"
                # self.push('puts "Disconnection by client request"\n')
                self.push('close\n')
                self.close()
                return             # necessary ?
            self.stack_count += 1
            self.push(c)
        else:
            for c in commands:
                self.stack_count += 1
                self.push(c+"\n")

        
    def getinput(self):
        print 'python->vmd> ',
        return sys.stdin.readline()


#
# Test code
#
if __name__ == '__main__':
    print """\
>>> [c = interactive('localhost')]
    done; then execute loop() to enter interactive mode:
>>> asyncore.loop()"""    
    c = interactive('localhost')

