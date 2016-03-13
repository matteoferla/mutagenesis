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


class foo:
    def __init__(self,text='text'):
        self._data='initialised value'

    @classmethod
    def bar(cls):
        x=cls.__new__(cls)
        x._data='gherkins'
        return x

    def __str__(self):
        return self._data


class foo2():
    def mymethod(self):
        self._data="method"

    def __init__(self):
        self.mymethod()

    def __str__(self):
        return self._data


def madness():
    s="ABC"
    for x in range(3):
        for y in range(3):
            if x>=y: #force x to be smaller
                continue
            else:
                print(s[0:x]+"X"+s[x+1:y]+"Y"+s[y+1:3])



if __name__ == "__main__":
    #print(foo2())
    madness()