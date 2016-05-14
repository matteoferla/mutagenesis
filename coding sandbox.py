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
        foo2.mymethod2(self)

    def __str__(self):
        return self._data

class bar:
    def __init__(self):
        self.a=1
        self.b=2
        self.c=3
        self.d=self.synonym(self.a)
        self.e=property(self.a, None, None, "I'm the 'x' property.")

    def synonym(self,original):
        def getx():
            return getattr(self, original)
        return property(getx, None, None, "I'm the 'x' property.")

        def mymethod2(self):
            self._data+="_method2"

class C:
    def __init__(self):
        self._x = None

    def getx(self):
        return self._x

    def setx(self, value):
        self._x = value

    def delx(self):
        del self._x

    x = property(getx,setx)


if __name__ == "__main__":
    #print(property(5, None, None, "I'm the 'x' property."))
    c=C()
    c.x=5
    print(c.x)
