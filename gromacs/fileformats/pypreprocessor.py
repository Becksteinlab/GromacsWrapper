#!/usr/bin/env python
# pypreprocessor.py

__author__ = 'Evan Plaice'
__version__ = '0.4.0'

import sys
import os
import traceback
import imp

class preprocessor:
    def __init__(self):
        # public variables
        self.defines = []
        self.input = sys.argv[0]
        self.output = ''
        self.removeMeta = False
        # private variables
        self.__linenum = 0
        self.__excludeblock = False
        self.__ifblock = False
        self.__ifcondition = ''
        self.__ifconditions = []
        self.__evalsquelch = True
        self.__outputBuffer = ''

    # the #define directive
    def define(self, define):
        self.defines.append(define)

    def search_defines(self, define):
        if define in self.defines:
            return True
        else:
            return False

    # the #ifdef directive
    def compare_defines_and_conditions(self, defines, conditions):
        # if defines and conditions lists have no intersecting values (ie. else = true)
        if not [val for val in defines if val in conditions]:
            return True
        else:
            return False

    # the #undef directive
    def undefine(self, define):
        # re-map the defines list excluding the define specified in the args
        self.defines[:] = [x for x in self.defines if x != define]

    # evaluate
    def lexer(self, line):
    # return values are (squelch, metadata)
        if self.__ifblock is False and self.__excludeblock is False:
            # squelch the preprocessor parse on the first
            # pass to prevent preprocessor infinite loop
            if 'pypreprocessor.parse()' in line:
                return True, True
            if line[:1] != '#':
                return False, False
        # handle #define directives
        if line[:7] == '#define':
            if len(line.split()) != 2:
                self.exit_error('#define')
            else:
                self.define(line.split()[1])
                return False, True
        # handle #undef directives
        if line[:6] == '#undef':
            if len(line.split()) != 2:
                self.exit_error('#undef')
            else:
                self.undefine(line.split()[1])
                return False, True
        # handle #endif directives
        if line[:6] == '#endif':
            if len(line.split()) != 1:
                self.exit_error('#endif')
            else:
                self.__ifblock = False
                self.__ifcondition = ''
                self.__ifconditions = []
                return False, True
        # handle #endexclude directives
        if line[:11] == '#endexclude':
            if len(line.split()) != 1:
                self.exit_error('#endexclude')
            else:
                self.__excludeblock = False
                return False, True
        # handle #exclude directives
        if line[:8] == '#exclude':
            if len(line.split()) != 1:
                self.exit_error('#exclude')
            else:
                self.__excludeblock = True
        # process the excludeblock
        if self.__excludeblock is True:
            return True, False
        # handle #ifdef directives
        if line[:6] == '#ifdef':
            if len(line.split()) != 2:
                self.exit_error('#ifdef')
            else:
                self.__ifblock = True
                self.__ifcondition = line.split()[1]
                self.__ifconditions.append(line.split()[1])
        # handle #else directives
        if line[:5] == '#else':
            if len(line.split()) != 1:
                self.exit_error('#else')
        # process the ifblock
        if self.__ifblock is True:
            # evaluate and process an #ifdef
            if line[:6] == '#ifdef':
                if self.search_defines(self.__ifcondition):
                    self.__evalsquelch = False
                else:
                    self.__evalsquelch = True
                return False, True
            # evaluate and process the #else
            elif line[:5] == '#else':
                if self.compare_defines_and_conditions(self.defines, self.__ifconditions):
                    self.__evalsquelch = False
                else:
                    self.__evalsquelch = True
                return False, True
            else:
                return self.__evalsquelch, False
        else:
            return False, False

    # error handling
    def exit_error(self, directive):
        print('File: "' + self.input + '", line ' + str(self.__linenum))
        print('SyntaxError: Invalid ' + directive + ' directive')
        sys.exit(1)
    def rewrite_traceback(self):
        trace = traceback.format_exc().splitlines()
        index = 0
        for line in trace:
            if index == (len(trace) - 2):
                print(line.replace("<string>", self.input))
            else:
                print(line)
            index += 1
    
    # parsing/processing
    def parse(self):
        # open the input file
        input_file = open(os.path.join(self.input),'r')
        try:
            # process the input file
            for line in input_file:
                self.__linenum += 1
                # to squelch or not to squelch
                squelch, metaData = self.lexer(line)
                # process and output
                if self.removeMeta is True: 
                    if metaData is True or squelch is True:
                        continue
                if squelch is True:
                    self.__outputBuffer += '#' + line
                    continue
                if squelch is False:
                    self.__outputBuffer += line
                    continue
        finally:
            input_file.close()
        self.post_process()

    # post-processor
    def post_process(self):
        try:
            # open file for output (no auto-run)
            if self.output != '':
                self.run = False
                output_file = open(self.output, 'w')           
            # open tmp file
            else:
                self.run = True
                self.output = 'tmp_' + os.path.basename(self.input)
                output_file = open(self.output, 'w')
            # write post-processed code to file
            output_file.write(self.__outputBuffer)
        finally:
            output_file.close()            
        # resolve postprocess stage depending on the mode
        if self.run == False:
            sys.exit(0)
        else:    
            # if this module is loaded as a library override the import
            if imp.lock_held() is True:
                    self.override_import()
            else:
                self.on_the_fly()
                # break execution so python doesn't
                # run the rest of the pre-processed code
                sys.exit(0)

    # postprocessor - override an import
    def override_import(self):
        try:
            moduleName = self.input.split('.')[0]
            tmpModuleName = self.output.split('.')[0]
            del sys.modules[moduleName]
            sys.modules[tmpModuleName] = __import__(tmpModuleName)
            sys.modules[moduleName] = __import__(tmpModuleName)
        except:
            self.rewrite_traceback()
        finally: 
            # remove tmp (.py & .pyc) files
            os.remove(self.output)
            os.remove(self.output + 'c')

    # postprocessor - on-the-fly execution
    def on_the_fly(self):
        try:
            exec(open(self.output,"rb").read())
        except:
            self.rewrite_traceback()
        finally:
            # remove tmp file
            os.remove(self.output)
     
pypreprocessor = preprocessor()
