import pyautogui as pg
import ctypes
#import pyperclip as pc
import time

def winshot():
    pg.keyDown("winleft")
    pg.hotkey("printscreen")
    pg.keyUp("winleft")
    time.sleep(0.5)
    pg.hotkey("right")

def get_region():
    if(fsize==False):
        while True:
           print("1st",pg.position())
           if ctypes.windll.user32.GetAsyncKeyState(0x01) == 0x8000:
                x0,y0 = pg.position()
                break
        print("end 1")

        time.sleep(1)

        while True:
            print("2nd",pg.position())
            if ctypes.windll.user32.GetAsyncKeyState(0x01) == 0x8000:
                x1,y1 = pg.position()
                break
        print("end 2")

        w = x1 - x0 
        h = y1 - y0

    elif(fsize==True):
        x0,y0=0,0
        w,h=pg.size()

    #print(x0,y0,w,h)
    return x0,y0,w,h
##
pageNum=300
restTime=6.0
fsize=False#FullSizeScreen

# ###
time.sleep(0)
here=get_region()
print("x0,y0,w,h",here)
print("start,after 5sec")
time.sleep(5)
for i in range(pageNum+1):
    #winshot()
    sc=pg.screenshot(region=here)
    pg.hotkey("right")
    filename="./img/img_%03d.png" %(i)
    #print(filename)
    sc.save(filename)
    print("page:",i)
    time.sleep(restTime)
print("end")







