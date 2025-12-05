class Functions(object):
    @staticmethod
    def format_print(lst: list):
        count = 0
        for num in lst:
            print("%3d" % num, end="  ")
            count += 1
            if (count % 8 == 0):
                print(" ")