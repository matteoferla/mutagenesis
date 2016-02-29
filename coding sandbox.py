__author__ = 'Matteo'
__doc__ = ''''''

N = "\n"
T = "\t"
# N="<br/>"

class test:
    def __init__(self,data):
        self.data=data
    def __getitem__(self,x):
        print(x)
        return self.data[x[0]][x[1]]


class outer:
    class inner:
        def __init__(self):
            a="ZORRO"
            print(a)
            outself=super()
            self.y=outself._data+" from "+a

        def __str__(self):
            return self.y

    def __init__(self):
        self._data="set in outer instance"
        print(self._data)
        #self._data=self.inner(self, "within")
        self._data=self.inner()
        print(self._data)


if __name__ == "__main__":
    print(outer()._data)