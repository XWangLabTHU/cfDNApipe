# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 09:55:10 2019

@author: zhang-wei
"""


from .Configure import Configure
import os
from hashlib import md5
from .cfDNA_utils import commonError, flatten, isAlphaOrDigit, rmEndString
import time
import sys
import subprocess
from multiprocessing import Pool
import pandas as pd

__metaclass__ = type


class StepBase:
    def __init__(self, stepNum=None, prevStep=None):
        if stepNum is not None:
            self.__stepID = stepNum
        elif (stepNum is None) and (prevStep is not None) and (prevStep is not True):
            self.__stepID = prevStep.getStepID() + 1
        else:
            self.__stepID = 1
        self.inputs = {}
        self.outputs = {}
        self.params = {}
        self.logpath = {}
        self.__isFinished = False
        self.__startTime = self.getCurTime()
        # attentionSteps is designed for command line program only
        self.attentionSteps = ["bismark", "bowtie2", "identifyAdapter"]

    # get stepID
    def getStepID(self,):
        return self.__stepID

    # set input parameters, all the input will be absolute Path
    def setInput(self, inputName, inputValue):
        if isinstance(inputName, list):
            if len(inputName) != len(inputValue):
                raise commonError("Number of input name and value not equal.")
            values = self.absolutePath(inputValue)
            for name, value in zip(inputName, values):
                self.inputs[name] = value
        else:
            self.inputs[inputName] = self.absolutePath(inputValue)

    # get input value by dict name (key)
    def getInput(self, inputName):
        return self.inputs[inputName]

    # get all input keys
    def getInputs(self,):
        return list(self.inputs.keys())

    # set output parameters
    def setOutput(self, outputName, outputValue):
        if isinstance(outputName, list):
            if len(outputName) != len(outputValue):
                raise commonError("Number of output key name and value not equal.")
            values = self.absolutePath(outputValue)
            for name, value in zip(outputName, values):
                self.outputs[name] = value
        else:
            self.outputs[outputName] = self.absolutePath(outputValue)

    # get output value by dict name (key)
    def getOutput(self, outputName):
        return self.outputs[outputName]

    # get all output keys
    def getOutputs(self,):
        return list(self.outputs.keys())

    # check input and output file path
    def checkFilePath(self, checkExist=True):
        self.checkOutputFilePath(checkExist)
        self.checkInputFilePath(checkExist)

    # check output file path
    def checkOutputFilePath(self, checkExist=True):
        for key in self.outputs.keys():
            outPaths = self.convertToList(self.outputs[key])
            self.checkFilePathList(outPaths, "output", key, checkExist)

    # check input file path
    def checkInputFilePath(self, checkExist=True):
        for key in self.inputs.keys():
            inPaths = self.convertToList(self.inputs[key])
            self.checkFilePathList(inPaths, "input", key, checkExist)

    # check path in a list
    def checkFilePathList(self, filePathList, iodict, key=None, checkExist=True):
        for filePath in filePathList:
            if filePath is None:
                raise commonError(
                    "File path of " + iodict + " " + key + " can not be None."
                )
            if checkExist:
                if not os.path.exists(filePath):
                    raise commonError(
                        "File path of "
                        + iodict
                        + " "
                        + key
                        + " not found: "
                        + filePath
                        + "."
                    )

    # get this Step folder name
    def getStepFolderName(self,):
        return "step_" + str(self.getStepID()).zfill(2) + "_" + self.__class__.__name__

    # get this step folder path
    def getStepFolderPath(self,):
        return Configure.getTmpPath(self.getStepFolderName())

    # create this step folder
    def stepFolderInit(self,):
        if not os.path.exists(self.getStepFolderPath()):
            os.mkdir(self.getStepFolderPath())

    # step init, include create step folder, log file and record file
    def stepInit(self, upstream=None):
        if upstream is not None:
            self.stepFolderInit()
            self.setPipeLogPath()
            self.setPipeRecPath()
        else:
            self.setLogPath(
                os.path.join(self.getOutput("outputdir"), self.getLogName())
            )
            self.setRecPath(
                os.path.join(self.getOutput("outputdir"), self.getRecName())
            )

        if os.path.exists(self.getLogPath()):
            finishFlag = self.checkFinish()
        else:
            finishFlag = False

        if not finishFlag:
            self.createLog(overwrite=True)
            self.createRec(overwrite=True)

        return finishFlag

    # get log file name
    def getLogName(self,):
        return self.__class__.__name__ + "." + self.getParaMD5code() + ".log"

    # set log path
    def setLogPath(self, path):
        self.logpath = path

    # set pipeline log path
    def setPipeLogPath(self,):
        self.logpath = os.path.join(self.getStepFolderPath(), self.getLogName())

    # get log file path
    def getLogPath(self,):
        return self.logpath

    # create log file
    def createLog(self, overwrite=False):
        if not os.path.exists(self.logpath):
            logs = {
                "Class_Name": {"value": "None", "Flag": False},
                "Start_Time": {"value": "None", "Flag": False},
                "End_Time": {"value": "None", "Flag": False},
                "Input_Files": {"value": "None", "Flag": False},
                "Output_Files": {"value": "None", "Flag": False},
                "Log_Files": {"value": "None", "Flag": False},
                "Record_Files": {"value": "None", "Flag": False},
                "Fi_Flag": {"value": "None", "Flag": False},
            }
            logs_df = pd.DataFrame.from_dict(logs, orient="index")
            logs_df.to_csv(self.getLogPath(), sep="\t")
        else:
            pass

    # load log
    def loadLog(self,):
        df = pd.read_csv(self.getLogPath(), sep="\t", index_col=0)
        dict_df = df.to_dict(orient="index")
        return dict_df

    # get log file value
    def getLogValue(self, key):
        dict_df = self.loadLog()
        request_value = dict_df[key]
        return request_value

    # save dict to log file
    def saveLog(self, dict_df):
        logs_df = pd.DataFrame.from_dict(dict_df, orient="index")
        logs_df.to_csv(self.getLogPath(), sep="\t")

    # get record file name
    def getRecName(self,):
        return self.__class__.__name__ + "." + self.getParaMD5code() + "-record.txt"

    # set record path
    def setRecPath(self, path):
        self.recpath = path

    # set pipeline record path
    def setPipeRecPath(self,):
        self.recpath = os.path.join(self.getStepFolderPath(), self.getRecName())

    # get log file path
    def getRecPath(self,):
        return self.recpath

    # create recort file
    def createRec(self, overwrite=False):
        if not os.path.exists(self.recpath):
            open(self.recpath, "a").close()
        else:
            if overwrite:
                open(self.recpath, "w").close()
            else:
                pass

    # write log
    def writeRec(self, mess):
        recfile = open(self.recpath, "a")
        recfile.writelines(mess + "\n")
        recfile.close()

    # check whether finished
    def checkFinish(self,):
        try:
            Fi_Flag = self.getLogValue("Fi_Flag")
        except KeyError:
            return False
        else:
            return Fi_Flag["Flag"]

    # set input and output parameters
    def setParam(self, paramName, paramValue):
        self.params[paramName] = paramValue

    # get input and output parameters
    def getParam(self, paramName):
        return self.params[paramName]

    # get all parameter keys
    def getParams(self,):
        return list(self.params.keys())

    # get time
    def getCurTime(self,):
        return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))

    # get file absolute path
    def absolutePath(self, pathOrPathList):
        if pathOrPathList is None:
            return None
        elif isinstance(pathOrPathList, list):
            return [os.path.abspath(s) for s in pathOrPathList]
        else:
            return os.path.abspath(pathOrPathList)

    # convert a odj to a list
    def convertToList(self, obj):
        if not isinstance(obj, list):
            return [obj]
        else:
            return obj

    # get single input file name
    def getFileNamePrefix(self, fileName):
        return os.path.split(fileName)[-1]

    # get paired input file name prefix
    def getMaxFileNamePrefix(self, file1, file2):
        file1 = self.getFileNamePrefix(file1)
        file2 = self.getFileNamePrefix(file2)
        len1 = len(file1)
        len2 = len(file2)
        for i in range(min(len1, len2)):
            if file1[i] != file2[i]:
                break
        if i == 0:
            return ""
        elif i == (min(len1, len2) - 1):
            tmp_str = file1[: i + 1]
        else:
            tmp_str = file1[:i]

        for k in reversed(range(len(tmp_str))):
            if isAlphaOrDigit(tmp_str[k]):
                final_name = tmp_str[: k + 1]
                final_name = rmEndString(final_name, [".pair"])
                return final_name
            else:
                k = k - 1

        raise commonError("File names must contain at least one alphabet or number.")

    # single prefix
    def getMaxFileNamePrefixV2(self, file):
        final_name = os.path.splitext(os.path.basename(file))[0]
        final_name = rmEndString(
            final_name,
            [
                "-sorted",  # bamsort suffix
                "_sorted",  # bamsort suffix
                "-rmdup",  # remove duplicates suffix
                ".fq.gz",
                ".fq",  # ***
                ".pair1.truncated",
                ".pair1.truncated.gz_bismark_bt2_pe",  # bisamrk WGBS paired suffix
                ".truncated.gz_bismark_bt2",  # bisamrk WGBS single suffix
                ".pair1.truncated.gz_unmapped_reads_1",
                ".lc_filter_R1",
                "_pe",
                "_se"
                ".pair1.truncated.gz_bismark_bt2_pe",  # bisamrk WGBS paired suffix
                "_pe.deduplicated",
                "-BQSR",
                "-RG",
                ".unfiltered.vcf",
                ".filtered.vcf",
                ".vcf",
                ".snp.raw",  # for BisSNP genotyper
                "_unmapped_reads_1",  # for bismark unmapped read1
                ".R1",  # for samtofast
                ".pair1.truncated.gz_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero"
                ".pair1.truncated.gz_bismark_bt2_pe.deduplicated",
                ".pair1.truncated.gz_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero",
                ".truncated.gz_bismark_bt2",  # bisamrk WGBS single suffix
                ".truncated.gz_bismark_bt2.deduplicated",
                ".truncated.gz_bismark_bt2.deduplicated.bedGraph.gz.bismark.zero",
            ],
        )
        return final_name

    # get file name and size
    def getFileNameAndSize(self, filePath, fileSize=True):
        namesizelist = []
        if not isinstance(filePath, list):
            filePath = [filePath]
        for filepath in filePath:
            filename = os.path.split(filepath)[-1]
            if fileSize:
                filesize = os.path.getsize(filepath)
            else:
                filesize = ""
            namesizelist.append(filename + str(filesize))
        return namesizelist

    # get md5 code, check input, output and parameters
    def getParaMD5code(self,):
        checklist = [""]
        checklist1 = []
        checklist2 = []
        checklist3 = []
        for key in self.inputs.keys():
            checklist1.extend(self.getFileNameAndSize(self.inputs[key]))
        checklist1.sort()
        checklist.extend(checklist1)
        for key in self.outputs.keys():
            checklist2.extend(
                self.getFileNameAndSize(self.outputs[key], fileSize=False)
            )
        checklist2.sort()
        checklist.extend(checklist2)
        keys = list(self.params.keys())
        keys.sort()
        for key in keys:
            checklist3.append(self.params[key])
        checklist.extend(str(checklist3))
        checklist = "".join(checklist)
        return md5(checklist.encode("utf8")).hexdigest()[0:8]

    # get configure value
    def getConfigVal(self, key):
        return Configure.getConfig(key)

    # command create, output a string
    def cmdCreate(self, cmdlist):
        if not isinstance(cmdlist, list):
            raise commonError("Parameter 'cmdlist' must be a list!")
        else:
            tmp_cmd = []
            for tmp_value in cmdlist:
                if isinstance(tmp_value, dict):
                    for k, v in tmp_value.items():
                        if isinstance(v, bool) and v:
                            tmp_cmd.append(k)
                        elif isinstance(v, bool) and (not v):
                            pass
                        else:
                            tmp_cmd.append(k)
                            tmp_cmd.append(v)
                elif isinstance(tmp_value, list):
                    for l in tmp_value:
                        tmp_cmd.append(l)
                else:
                    tmp_cmd.append(tmp_value)

        cmd = " ".join([str(x) for x in tmp_cmd])
        return cmd

    # run the command line, this function is designed for sngle threads
    def run(self, cmds, force=False):
        self.writeRec(
            "#############################################################################################"
        )
        # read log file and get any executed cmd
        logs_dict = self.loadLog()
        fi_cmd = []
        for tmp_key in logs_dict:
            if tmp_key.startswith("cmd_"):
                if logs_dict[tmp_key]["Flag"]:
                    fi_cmd.append(logs_dict[tmp_key]["value"])
        if isinstance(cmds, list):  # cmd is a list
            for idx, cmd in enumerate(cmds):
                # whether the cmd is finished
                if (cmd in fi_cmd) and not force:
                    mess = "#" * 20
                    print(mess)
                    print("The following command had been completed!")
                    print(cmd)
                    print(mess)
                    continue
                self.writeRec("Cmd: {}".format(cmd))
                print("Now, running command: {}".format(cmd))
                proc = subprocess.Popen(
                    cmd,
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                )
                while True:
                    nextline = proc.stdout.readline()
                    if (nextline == "") and (proc.poll() is not None):
                        break
                    sys.stdout.write(nextline)
                    self.writeRec(nextline.strip())
                    sys.stdout.flush()

                output, error = proc.communicate()
                exitCode = proc.returncode

                # catch error
                if exitCode == 0:  # successfully finished
                    if self.__class__.__name__ in self.attentionSteps:
                        cmd_key = "cmd_" + md5(cmd.encode("utf-8")).hexdigest()
                        logs_dict[cmd_key] = {"value": cmd, "Flag": True}
                        self.saveLog(logs_dict)
                else:
                    self.writeRec("ERROR CMD: {}".format(cmd))
                    self.writeRec("ERROR CODE: {}".format(exitCode))
                    self.writeRec("ERROR MESSAGE: {}".format(error))
                    self.writeRec(
                        "#############################################################################################"
                    )
                    if self.__class__.__name__ in self.attentionSteps:
                        cmd_key = "cmd_" + md5(cmd.encode("utf-8")).hexdigest()
                        logs_dict[cmd_key] = {"value": cmd, "Flag": False}
                        self.saveLog(logs_dict)

                    raise commonError("**********CMD running error**********")

        else:  # cmd is a string
            self.writeRec("Cmd: {}".format(cmds))
            print("Now, running command: {}".format(cmds))
            proc = subprocess.Popen(
                cmds,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
            )
            while True:
                nextline = proc.stdout.readline()
                if (nextline == "") and (proc.poll() is not None):
                    break
                sys.stdout.write(nextline)
                self.writeRec(nextline.strip("\n"))
                sys.stdout.flush()

            output, error = proc.communicate()
            exitCode = proc.returncode

            # catch error
            if exitCode == 0:  # successfully finished
                if self.__class__.__name__ in self.attentionSteps:
                    cmd_key = "cmd_" + md5(cmds.encode("utf-8")).hexdigest()
                    logs_dict[cmd_key] = {"value": cmds, "Flag": True}
                    self.saveLog(logs_dict)
            else:
                self.writeRec("ERROR CMD: {}".format(cmds))
                self.writeRec("ERROR CODE: {}".format(exitCode))
                self.writeRec("ERROR MESSAGE: {}".format(error))
                self.writeRec(
                    "#############################################################################################"
                )
                if self.__class__.__name__ in self.attentionSteps:
                    cmd_key = "cmd_" + md5(cmds.encode("utf-8")).hexdigest()
                    logs_dict[cmd_key] = {"value": cmds, "Flag": False}
                    self.saveLog(logs_dict)

                raise commonError("**********CMD running error**********")

    # step information record, every step will use this method
    def stepInfoRec(self, cmds, finishFlag):
        logs_dict = self.loadLog()
        self.setParam("cmd", list(flatten(cmds)))
        if finishFlag:
            mess = "*" * 60
            print(mess)
            print("{:^60s}".format(self.__class__.__name__ + " has been completed!"))
            print(mess)
        else:
            logs_dict["Class_Name"]["value"] = self.__class__.__name__
            logs_dict["Start_Time"]["value"] = self.__startTime
            logs_dict["End_Time"]["value"] = self.getCurTime()
            logs_dict["Input_Files"]["value"] = self.inputs.values()
            logs_dict["Output_Files"]["value"] = self.outputs.values()
            logs_dict["Log_Files"]["value"] = self.getLogPath()
            logs_dict["Record_Files"]["value"] = self.getRecPath()
            logs_dict["Fi_Flag"]["Flag"] = True

            # cmd information
            for idx, cmd in enumerate(self.getParam("cmd")):
                tmp_key = "CMD_No." + md5(cmd.encode("utf-8")).hexdigest()
                logs_dict[tmp_key] = {"value": cmd, "Flag": True}

            self.writeRec("Record finished!")
            self.saveLog(logs_dict)

    # single function run, designed for multicore
    def funRun(self, args):
        try:
            fun = args[0]
            param = args[1]
            results = fun(*param)
            flag = True
        except Exception as e:
            results = e
            flag = False

        return results, flag

    # single command run, designed for multicore
    def cmdRun(self, cmd):
        print(cmd)
        proc = subprocess.run(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        # get cmd info
        mess = proc.stderr + proc.stdout
        if proc.returncode == 0:
            return mess, True
        else:
            mess1 = (
                "\n\n\nAn Error Occured During The Following Command Line Executing.\n"
            )
            mess2 = (
                "\n         Please Stop The Program To Check The Error.         \n\n\n"
            )
            mess = """^^^{}^^^\n{}\n^^^{}^^^""".format(mess1, cmd, mess2)
            return mess, False

    # time counter for multiRun
    def track_job(self, job, time_start, update_interval=3, print_interval=300):
        """
        from stackoverflow
        job: multiprocessing job
        """
        minute_count = 0
        while job._number_left > 0:
            delta_minutes = (time.time() - time_start) // print_interval
            if delta_minutes and (delta_minutes != minute_count):
                print(
                    "{0} minutes has passed since the last print".format(
                        delta_minutes * 5
                    )
                )
                print("Tasks remaining = {0}".format(job._number_left * job._chunksize))
                minute_count = delta_minutes

            time.sleep(update_interval)

    # multiCore
    def multiRun(self, args, func=None, nCore=1):
        """
        This function is designed for multiCore.

        multiRun(func, args, type, nCore=1)
        {P}arameters:
            args: Parameters for multirun, [[1, 2, 3], [2, 3, 4]] or ["cmd1", "cmd2", "cmd3"].
            func: function name for multicore, None mean cmd mode.
            nCore: int, how many cores to use.
        """
        # time start
        time_start = time.time()
        # load log file
        logs_dict = self.loadLog()

        p = Pool(nCore)
        print("Start multicore running, master process number: {}".format(nCore))
        print(
            "Note: some command line verbose may be blocked, the program will record them in record file."
        )

        if func is None:
            # in this mode, success means that the output will be stdout, failed means that the output will be False
            if self.__class__.__name__ in self.attentionSteps:
                fi_cmd = []
                for tmp_key in logs_dict:
                    if tmp_key.startswith("cmd_"):
                        if logs_dict[tmp_key]["Flag"]:
                            fi_cmd.append(logs_dict[tmp_key]["value"])
                # print finished command line
                print("#" * 20)
                for idx, cmd in enumerate(args):  # whether the cmd is finished
                    if cmd in fi_cmd:
                        print(cmd + " had been completed!")
                        args.remove(cmd)
                        continue
                print("#" * 20)

            results = p.map_async(self.cmdRun, args)
        else:
            # in this mode, success means that the output will be function output
            reshaped_args = [[[func, x]] for x in args]
            results = p.starmap_async(self.funRun, reshaped_args)

        # print mess
        print("Subprocesses Start running......")
        print("Waiting for all subprocesses done...")

        self.track_job(
            job=results, time_start=time_start, update_interval=3, print_interval=300
        )

        p.close()
        p.join()
        print("All subprocesses done.")

        if func is None:  # check output for cmd
            messages = [x[0] for x in results.get()]
            flag = [x[1] for x in results.get()]
            mess = "\n\n\n".join([str(x) for x in messages])

            for i, j in zip(args, flag):
                cmd_key = "cmd_" + md5(i.encode("utf-8")).hexdigest()
                logs_dict[cmd_key] = {"value": i, "Flag": j}
                self.saveLog(logs_dict)

            if (all(flag)) and results.successful():
                return mess, True
            else:
                print(mess)
                raise commonError("Error occured in multi-core running!")
        else:  # check output for function
            output = [x[0] for x in results.get()]
            flag = [x[1] for x in results.get()]
            if (all(flag)) and results.successful():
                return output, True
            else:
                print(str(output))
                raise commonError("Error occured in multi-core running!")
