# VSCode intro

VSCode is a powerful application for programming. Fundamentally it is a text editor, but it additionally has all sorts of extra plugins that allow it to function as a powerful development environment.

Download VSCode here for your system: https://code.visualstudio.com/download

Importantly, VSCode can connect to Sherlock and allow you to edit your Sherlock files without the use of vim/emacs/nano/whatever. When connecting to a remote server, VSCode downloads the files locally for editing, and then reuploads them back to the remote server.

This is how it looks for me when I edit code:
<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/75d9bdac-7956-485e-bebf-db6c0ca81a28" />

What is really going on here is VSCode is accessing sherlock and downloading my files (specifically the ones I open to edit, not everything) on to my local computer. It then opens the file and allows me to edit using the VSCode interface. When I then save any changes that have been made, VSCode automatically uploads the changed file back to sherlock, where it can be run.


### Setup SFTP

In order to connect to Sherlock, we will need to first create a directory somewhere on our local computer for VSCode to download files into (ie. ~/Documents/sherlock/). Next we will need to go to the Extensions tab (Fifth tab on the extreme left side of the window; Cmd/Ctrl+Shift+X) and download and install the SFTP plugin by publisher:"Natizyskunk". 

The extension page has instructions on how to connect to a remote server which I will summarize here.

First we need to open the local directory in VSCode so that it knows we will want to download files specifically into that location. Open the Explorer tab (First tab on extreme left side; Cmd/Ctrl+Shift+E) and hit 'Open Folder' and select the new directory.

Next go to the SFTP tab (Should be the last tab on the extreme left; File and cloud icon). We need to set up a config file that tells VSCode what to connect to. Hit Cmd/Ctrl+Shift+P and type in `SFTP: Config`/ which should open a file called `sftp.json`.

Here is the config file I use for Sherlock. This config file will give me two separate tabs in the SFTP window. One will connect me to `$SCRATCH` on Sherlock and the other will connect me to `$HOME`:
```
[
{
    "name": "sherlock_home",
    "context": "sherlock_home",
    "host": "dtn.sherlock.stanford.edu",
    "protocol": "sftp",
    "port": 22,
    "username": "jahemker",
    "remotePath": "/home/users/jahemker/",
    "uploadOnSave": true
},
{
    "name": "sherlock_scratch",
    "context": "sherlock_scratch",
    "host": "dtn.sherlock.stanford.edu",
    "protocol": "sftp",
    "port": 22,
    "username": "jahemker",
    "remotePath": "/scratch/users/jahemker/",
    "uploadOnSave": true
}
]
```
`name` and `context` define what it is we are connecting to. 

`host` is the address connection essentially. We want to use `dtn.sherlock.stanford.edu` instead of `login.sherlock.stanford.edu` because all we are doing is asking to transfer files. The `dtn` nodes will ask for you password, but they won't ask for two-factor authentication.

`protocol` and `port` should just be `sftp` and a number (22 is fine).

`username` is your username to log into Sherlock.

`remotePath` is where the sftp connection will place you in the server. You can always go deeper in the tree, but you can't zoom out, so I would recommend just having this as the base-most directory, which is `$HOME` and `$SCRATCH` on Sherlock.

`uploadOnSave` can be true or false. I prefer it to be true. This means that everytime you save a file you are locally editing it will automatically upload the new version back to Sherlock. I think if it's false you just have to manually trigger the upload back to Sherlock. I believe you would do this by hitting Cmd/Ctrl+Shift+P and using the command `SFTP: Upload Active File`.

There is a `password` field that you can use and put your Sherlock password into so that you don't have to type it in each session, but it doesn't encrypt it at all so I would not recommend doing that.

All configuration settings can be found here: https://github.com/Natizyskunk/vscode-sftp/wiki/Common-Configuration

---

You should now be able to see a couple drop down menus named `sherlock_home` and `sherlock_scratch` within the SFTP tab. You can basically treat this as a mini file explorer for Sherlock. You should be able to open files, edit and save them, and see the changes on Sherlock.

<img width="301" height="331" alt="image" src="https://github.com/user-attachments/assets/f726e251-dea5-4597-af5b-96090b3350c8" />


**The SFTP explorer does not auto refresh**, so if you are making new directories or files on Sherlock, you will need to hit the little refresh button to see them in VSCode.

Sometimes when VSCode is restarted, the SFTP extension tab will not be present and it will ask you to choose a workspace directory. You just need to point it in the direction of the original directory you set up for sherlock on your local computer (in this token example it was ~/Documents/sherlock/).

