import logging

class LoggerForVaspTask:

    @staticmethod
    def log_output(**kwargs):
        level = kwargs.pop("level", None)
        log_file_path = kwargs.pop("log_file_path", None)
        FORMAT = kwargs.pop("FORMAT", None)
        DATEFMT = kwargs.pop("DATEFMT", None)
        msg = kwargs.pop("msg", None)

        assert ( type(level) == int )    
        assert ( type(log_file_path) == str )
        assert ( type(msg) == str )

        if level is None:
            level = logging.DEBUG
        if FORMAT is None:
            FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        if DATEFMT is None:
            DATEFMT = "%Y-%m-%d %H:%M:%S"

        logger = logging.getLogger("VASP_Task_logger")
        logger.setLevel(level)
        logger.propagate = False

        formatter = logging.Formatter(FORMAT, DATEFMT)
        
        sh = logging.FileHandler(log_file_path)
        sh.setLevel(level)
        sh.setFormatter(formatter)

        logger.addHandler(sh)

       
        if logger.level == logging.DEBUG:
            logger.debug(msg)
        elif logger.level == logging.INFO:
            logger.info(msg)
        elif logger.level == logging.WARNING:
            logger.warning(msg)
        elif logger.level == logging.ERROR:
            logger.error(msg)
        elif logger.level == logging.CRITICAL:
            logger.critical(msg)

        # 关闭 log 文件
        sh.close()
        logger.removeHandler(sh)