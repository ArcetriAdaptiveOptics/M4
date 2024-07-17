'''
HOW TO USE IT::

    from m4.ground.timestamp import Timestamp
    now = Timestamp.now()
    today = Timestamp.today()
    or
    ts = Timestamp()
    now = ts.asNowString()
    today = ts.asTodayString()
'''
import datetime


class Timestamp():
    ''' Class for tracking numbers generation
    '''

    def __init__(self, tn=None):
        if tn is None:
            self._now = datetime.datetime.now()
        else:
            self._now = datetime.datetime.strptime(
                tn, '%Y%m%d_%H%M%S')

    def asNowString(self):
        '''
        Returns
        -------
        now: string
            yearMonthDay_hoursMinutesSeconds
        '''
        return self._now.strftime("%Y%m%d_%H%M%S")

    def asTodayString(self):
        '''
        Returns
        -------
        today: string
            year month day
        '''
        return self._now.strftime("%Y%m%d")

    @staticmethod
    def now():
        '''
        Returns
        -------
        now: string
            yearMonthDay_hoursMinutesSeconds
        '''
        return Timestamp().asNowString()

    @staticmethod
    def today():
        '''
        Returns
        -------
        today: string
            year month day
        '''
        return Timestamp().asTodayString()

    @staticmethod
    def fromString(tn):
        '''
        Returns
        -------
        today: Timestamp object
        '''
        return Timestamp(tn)

    def __str__(self):
        return self.asNowString()

    @property
    def datetime(self):
        return self._now
