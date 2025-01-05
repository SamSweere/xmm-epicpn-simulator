import os
import subprocess
from pathlib import Path

import heasoftpy as hsp
from loguru import logger

_default_args = {
    "noprompt": True,
    "stderr": True,
    "chatter": 0,
}


class HSPTask(hsp.HSPTask):
    def __call__(self, args=None, **kwargs) -> hsp.HSPResult:
        """Call the task.

        There are several ways to call HSPTask:
        1: HSPTask(): required parameters will be queried as usually done with heasoft tasks.
        2: HSPTask(previousHSPTask): initialize from a previously defined HSPTask.
        3: HSPTask(argsDict): argsDict is a dict of input parameters.
        4: HSPTask(foo=bar, x=0): parameters are passed in the kwargs dict.


        Parameters
        ----------
        args: HSPTask or dict
            Task parameters as another HSPTask or dict

        Keywords
        --------
        individual task parameters given as: paramter=value.

        Additionally, the user may pass other arguments that are common between
        all tasks. These include:

        verbose: This can take several values. In all cases, the text printed by the
            task is captured, and returned in HSPResult.stdout/stderr. Addionally:
            - 0 (also False or 'no'): Just return the text, no progress prining.
            - 1 (also True or 'yes'): In addition to capturing and returning the text,
                task text will printed into the screen as the task runs.
            - 2: Similar to 1, but also prints the text to a log file.
            - 20: In addition to capturing and returning the text, log it to a file,
                but not to the screen.
                In both cases of 2 and 20, the default log file name is {taskname}.log.
                A logfile parameter can be passed to the task to override the file name.

        noprompt: Typically, HSPTask would check the input parameters and
            queries any missing ones. Some tasks (e.g. pipelines) can run by using
            default values. Setting `noprompt=True`, disables checking and quering
            the parameters. Default is False.

        stderr: If True, make stderr separate from stdout. The default
            is False, so stderr is written to stdout.

        Returns:
            HSPResult instance.
            e.g. HSPResult(ret_code, std_out, std_err, params, custom_dict)
        """

        # assemble the user input, if any, into a dict
        if args is None:
            user_pars = {}
        elif isinstance(args, dict):
            user_pars = dict(args)
        elif isinstance(args, HSPTask):
            user_pars = dict(args.params)
        else:
            raise hsp.HSPTaskException("Unrecognized input in initializing HSPTask")

        # also any parameters in self.params from a previous call
        # or entered by hand
        user_pars.update(self.params)

        # add all keywords if present
        # any commandLine arguments in sys.argv should have already been processed into kwargs
        user_pars.update(kwargs)

        # ----------------------------- #
        # handle common task parameters #
        # do we have an explicit stderr
        stderr = user_pars.get("stderr", False)
        if not isinstance(stderr, bool):
            stderr = (
                (isinstance(stderr, str) and stderr.strip().lower() in ["y", "yes", "true"])
                or isinstance(stderr, int)
                and stderr > 0
            )
        self.stderr = stderr

        # noprompt?
        noprompt = user_pars.get("noprompt", False)
        if "noprompt" in self.par_names:
            # in case noprompt is task parameter, we look for py_noprompt
            noprompt = user_pars.get("py_noprompt", False)
        if not isinstance(stderr, bool):
            noprompt = (
                (isinstance(noprompt, str) and noprompt.strip().lower() in ["y", "yes", "true"])
                or isinstance(noprompt, int)
                and noprompt > 0
            )
        self._noprompt = noprompt

        # verbose?
        verbose = user_pars.get("verbose", 0)
        if "verbose" in self.par_names:
            # in case verbose is task parameter, we look for py_verbose
            verbose = user_pars.get("py_verbose", 0)

        if isinstance(verbose, bool):
            verbose = 1 if verbose else 0
        elif isinstance(verbose, str):
            if verbose.strip().lower() in ["y", "yes", "true"]:
                verbose = 1
            elif verbose.strip().lower() in ["n", "no", "false"]:
                verbose = 0
            else:
                try:
                    verbose = int(verbose)
                except ValueError:
                    verbose = 1
        if not isinstance(verbose, int):
            raise hsp.HSPTaskException("confusing verbose value. Allowed types are: bool, str or int")
        self._verbose = verbose

        # now check the user input against expectations, and query if incomplete
        usr_params = self.build_params(user_pars)

        # create a dict for all model parameters
        params = {p: getattr(self, p).value for p in self.par_names}
        self.params = usr_params if self._noprompt else params

        # disable prompt: https://heasarc.gsfc.nasa.gov/lheasoft/scripting.html
        os.environ["HEADASNOQUERY"] = ""
        os.environ["HEADASPROMPT"] = "/dev/null"

        # write new params to the user .par file
        # do this before calling in case the task also updates the .par file
        usr_pfile = HSPTask.find_pfile(self.taskname, return_user=True)
        self.write_pfile(usr_pfile)

        # now call the task #
        result = self.exec_task()

        # ensure we are returning the correct type
        if not isinstance(result, hsp.HSPResult):
            raise hsp.HSPTaskException(f"Returned result type {type(result)} is not HSPResult")

        # re-read the pfile in case it has been modified by the task
        # update only the the values in the HSPTask instance, not
        #  result.params that will be returned to the user
        if os.path.exists(usr_pfile):
            params_after = HSPTask.read_pfile(usr_pfile)
            for ipar, par_name in enumerate(self.par_names):
                setattr(self, par_name, params_after[ipar].value)
            # result.params.update(self.params)

        return result

    def exec_task(self):
        """Run the Heasoft task, but in our own custom way.

        Returns
        -------
        HSPResult instance.
            e.g. HSPResult(ret_code, std_out, std_err, params, custom_dict)

        """

        # Get the task parameters
        usr_params = self.params

        # construct a parameter list
        for par in usr_params:
            if isinstance(usr_params[par], bool):
                usr_params[par] = "yes" if usr_params[par] else "no"
            if usr_params[par] is None:
                usr_params[par] = "NONE"

            # '$( )' ensures empty string are passed correctly with subprocess
            if isinstance(usr_params[par], str) and usr_params[par] == "":
                usr_params[par] = "$( )"
        cmd_params = [f"{par}={val}" for par, val in usr_params.items()]

        # the task executable
        exec_cmd = os.path.join(os.environ["HEADAS"], f"bin/{self.taskname}")

        if os.path.exists(exec_cmd):
            exec_cmd = [exec_cmd]
        elif os.path.exists(exec_cmd + ".py"):
            exec_cmd = ["python", exec_cmd + ".py"]
        else:
            raise hsp.HSPTaskException(f"There is no task file {exec_cmd} to run")

        cmd_list = exec_cmd + cmd_params

        proc = subprocess.run(
            args=cmd_list,
            env=os.environ.copy(),
            capture_output=True,
            text=True,
            timeout=3600,
        )

        proc_out, proc_err = proc.stdout, proc.stderr
        if isinstance(proc_out, bytes):
            proc_out = proc_out.decode("ISO-8859-15")
        if isinstance(proc_err, bytes):
            proc_err = proc_err.decode("ISO-8859-15")

        return hsp.HSPResult(proc.returncode, proc_out, proc_err, usr_params)


def _check_result(task_name: str, hsp_result: hsp.HSPResult) -> None:
    if hsp_result.returncode != 0 or hsp_result.stderr:
        logger.error(f"{task_name} failed with return code {hsp_result.returncode} and:\n{hsp_result}")
        raise hsp.HSPTaskException(f"{task_name} failed with return code {hsp_result.returncode} and:\n{hsp_result}")

    logger.debug(f"{task_name} output:\n{hsp_result}")


def ftmerge(args=None, **kwargs) -> Path:
    if args is None:
        args = _default_args
    else:
        args.update(_default_args)

    with hsp.utils.local_pfiles_context():
        ftmerge_task = HSPTask(name="ftmerge")
        result = ftmerge_task(args, **kwargs)

    _check_result("ftmerge", result)

    return Path(result.params["outfile"])


def ftcopy(args=None, **kwargs) -> Path:
    if args is None:
        args = _default_args
    else:
        args.update(_default_args)

    with hsp.utils.local_pfiles_context():
        ftcopy_task = HSPTask(name="ftcopy")
        result = ftcopy_task(args, **kwargs)

    _check_result("ftcopy", result)

    return Path(result.params["outfile"])


def fthedit(args=None, **kwargs) -> None:
    if args is None:
        args = _default_args
    else:
        args.update(_default_args)

    with hsp.utils.local_pfiles_context():
        fthedit_task = HSPTask(name="fthedit")
        result = fthedit_task(args, **kwargs)

    _check_result("fthedit", result)


def ftedit(args=None, **kwargs) -> None:
    if args is None:
        args = _default_args
    else:
        args.update(_default_args)

    with hsp.utils.local_pfiles_context():
        ftedit_task = HSPTask(name="ftedit")
        result = ftedit_task(args, **kwargs)

    _check_result("ftedit", result)


def ftimgcalc(args=None, **kwargs) -> Path:
    if args is None:
        args = _default_args
    else:
        args.update(_default_args)

    with hsp.utils.local_pfiles_context():
        ftimgcalc_task = HSPTask(name="ftimgcalc")
        result = ftimgcalc_task(args, **kwargs)

    _check_result("ftimgcalc", result)

    return Path(result.params["outfile"])


def fimgbin(args=None, **kwargs) -> Path:
    if args is None:
        args = _default_args
    else:
        args.update(_default_args)

    with hsp.utils.local_pfiles_context():
        fimgbin_task = HSPTask(name="fimgbin")
        result = fimgbin_task(args, **kwargs)

    _check_result("fimgbin", result)

    return Path(result.params["outfile"])


def fimgmerge(args=None, **kwargs) -> Path:
    if args is None:
        args = _default_args
    else:
        args.update(_default_args)

    with hsp.utils.local_pfiles_context():
        fimgmerge_task = HSPTask(name="fimgmerge")
        result = fimgmerge_task(args, **kwargs)

    _check_result("fimgmerge", result)

    return Path(result.params["outfile"])
