var cfacdbOverlay = {
    
    win: null,
    
    openDialogCB: function(event)
    {
        if (!this.win || this.win.closed) {
            this.win = window.open('chrome://cfacdb/content/cfacdb.xul',        
                        'CFACDB', 'chrome,resizable=yes');
        } else {
            this.win.focus();
        }
    }
};
