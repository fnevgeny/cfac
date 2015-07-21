Components.utils.import("resource://gre/modules/AddonManager.jsm");

var cfacdb = {
    ini_id: -1,
    fin_id: -1,
    
    dnele: 0,
    
    dsources: null,
    
    logger: function(level, msg)
    {
        var prefix;

        switch (level) {
        case "err":
            prefix = "E";
            break;
        case "wrn":
            prefix = "W";
            break;
        case "inf":
            prefix = "I";
            break;
        case "dbg":
            if (!this.logDebug) {
                return;
            }
            prefix = "D";
            break;
        default:
            prefix = "?";
            break;
        }

        this.alert(prefix + ": " + msg);
    },
    
    setClassParams: function(class_name, value)
    {
        var es = document.getElementsByClassName(class_name);
        for (var i = 0; i < es.length; i++) {
            es[i].textContent = value;
        }
    },
        
    refresh: function(dsources)
    {
        if (!dsources) {
            return false;
        }
        
        var e;
        e = document.getElementById("levels-ini");
        e.datasources = dsources;
        
        e = document.getElementById("levels-fin");
        e.datasources = dsources;
        
        e = document.getElementById("rtransitions");
        e.datasources = dsources;
        
        e = document.getElementById("aitransitions");
        e.datasources = dsources;
        
        e = document.getElementById("ctransitions");
        e.datasources = dsources;
        
        e = document.getElementById("dE");
        e.datasources = dsources;
        
        e = document.getElementById("species");
        e.datasources = dsources;
        
        if (e.view && e.view.rowCount) {
            this.dsources = dsources;
            e.view.selection.select(0);
            return true;
        } else {
            return false;
        }
    },
    
    open: function(db_path)
    {
        var Cc = Components.classes;
        var Ci = Components.interfaces;
        
        var osString = Cc["@mozilla.org/xre/app-info;1"]
                         .getService(Ci.nsIXULRuntime).OS;
                   
        var dsources;
        if (osString == "WINNT") {
            dsources = "file:///" + db_path.replace(/\\/g, "/");
        } else {
            dsources = "file://" + db_path;
        }
        
        if (!this.refresh(dsources)) {
            this.alert("Invalid or inaccessible database!");
        }
    },

    openCB: function()
    {
        var Cc = Components.classes;
        var Ci = Components.interfaces;

        const nsIFilePicker = Ci.nsIFilePicker;

        var fp;
        if (!this.fpicker) {
            fp = Cc["@mozilla.org/filepicker;1"].createInstance(nsIFilePicker);
            fp.init(window, "Open File", nsIFilePicker.modeOpen);
            fp.appendFilter("cFAC SQLite databases", "*.sqlite; *.db");
            fp.appendFilters(nsIFilePicker.filterAll);
            this.fpicker = fp;
        } else {
            fp = this.fpicker;
        }

        fp.displayDirectory = this.lastDir;

        var rv = fp.show();

        if (rv == nsIFilePicker.returnOK) {
            this.setBusyCursor(true);
            
            var res;
            
            var file = fp.file;
            this.lastDir = file.parent;
            
            this.open(file.path);
            
            this.setBusyCursor(false);
            
            return res;
        } else {
            return false;
        }
        
    },
    
    sessionSelectCB: function(e)
    {
        var tree = e.target;
       
        var sid = 0;
        var nele_min = 0, nele_max = 1;
        
        if (tree.currentIndex >= 0) {
            var selection = tree.view.selection;
            var cell_text;

            cell_text = tree.view.getCellText(tree.currentIndex,  
                                tree.columns.getNamedColumn("sid"));

            if (cell_text) {
                sid = parseInt(cell_text);
            } else {
                sid = 0;
            }

            cell_text = tree.view.getCellText(tree.currentIndex,  
                                tree.columns.getNamedColumn("nele_min"));
            if (cell_text) {
                nele_min = parseInt(cell_text);
            } else {
                nele_min = 0;
            }
            cell_text = tree.view.getCellText(tree.currentIndex,  
                                tree.columns.getNamedColumn("nele_max"));
            if (cell_text) {
                nele_max = parseInt(cell_text);
            } else {
                nele_max = 1;
            }

            var nele_e = document.getElementById("nele");
            nele_e.min = nele_min;
            nele_e.max = nele_max;
        }
        
        this.setClassParams("sid", sid);
        
        var el;
        el = document.getElementById("levels-ini");
        el.view.selection.clearSelection();
        el.builder.rebuild();
        el = document.getElementById("levels-fin");
        el.view.selection.clearSelection();
        el.builder.rebuild();
    },
    
    levelSelectCB: function(e)
    {
        var tree = e.target;
       
        var id;
        if (tree.currentIndex < 0) {
            id = -1;
        } else {
            var selection = tree.view.selection;
            var cellText = tree.view.getCellText(tree.currentIndex,  
                                tree.columns.getPrimaryColumn());
            if (cellText) {
                id = parseInt(cellText);
            } else {
                id = -1;
            }
        }
        
        var class_name;
        
        if (tree.id == "levels-ini") {
            this.ini_id = id;
            class_name = "ini_id";
        } else {
            this.fin_id = id;
            class_name = "fin_id";
        }
        
        this.setBusyCursor(true);
        
        this.setClassParams(class_name, id);
        
        var el;
        var els = [];
        
        el = document.getElementById("dE");
        el.builder.rebuild();
        
        if (this.dnele == 0) {
            el = document.getElementById("rtransitions");
        } else {
            el = document.getElementById("aitransitions");
        }
        els.push(el);
        el = document.getElementById("ctransitions");
        els.push(el);
        
        for (var i in els) {
            el = els[i];
            if (this.fin_id < 0 || this.ini_id < 0) {
                el.datasources = "rdf:null";
            } else {
                if (el.datasources != this.dsources) {
                    el.datasources = this.dsources;
                } else {
                    // explicitly call rebuild() if the datasources haven't
                    // changed
                    el.builder.rebuild();
                }
            }
        }
        
        this.setBusyCursor(false);
    },
    
    neleSelectCB: function(e)
    {
        var el = e.target;
        var nele = parseInt(el.value);
        
        this.setClassParams("nele", nele);
        
        document.getElementById("levels-ini").builder.rebuild();
        document.getElementById("levels-fin").builder.rebuild();
    },
    
    deltaNeleCB: function(e)
    {
        var el = e.target;
        this.dnele = parseInt(el.value);
        this.setClassParams("dnele", this.dnele);
        
        document.getElementById("transitions-deck").selectedIndex = this.dnele;
        
        if (this.dsources) {
            var el = document.getElementById("levels-fin");
            if (el.selection) {
                el.selection.clearSelection();
            }
            el.builder.rebuild();
        }
    },
    
    refreshCB: function()
    {
        this.refresh(this.dsources);
    },
    
    closeCB: function()
    {
        window.close();
    },

    alert: function(text)
    {
        var Cc = Components.classes;
        var Ci = Components.interfaces;

        var prompts = Cc["@mozilla.org/embedcomp/prompt-service;1"]
                        .getService(Ci.nsIPromptService);
        return prompts.alert(null, "CFACDB - Alert", text);
    },

    aboutCB: function()
    {
        var httpBASE = this.httpBASE;
        AddonManager.getAddonByID("cfacdb@plasma-gate.weizmann.ac.il",
            function(aAddon) {
                var str = aAddon.name + " version: " + aAddon.version;
                str += "\nCreated by: " + aAddon.creator.name;
                cfacdb.alert(str);
            });
    },
    
    helpCB: function()
    {
        var url = "https://www-amdis.iaea.org/CFACDB/";
        var gBrowser = window.opener.gBrowser;

        if (gBrowser) {
            gBrowser.selectedTab = gBrowser.addTab(url);  
        }
    },
    
    aboutCB: function()
    {
        var httpBASE = this.httpBASE;
        AddonManager.getAddonByID("cfacdb@plasma-gate.weizmann.ac.il",
            function(aAddon) {
                var str = aAddon.name + " version: " + aAddon.version;
                str += "\nCreated by: " + aAddon.creator.name;
                // str += "\nREST Base: " + httpBASE;
                cfacdb.alert(str);
            });
    },
    
    setBusyCursor: function(onoff)
    {
        window.setCursor(onoff ? "wait" : "auto");
    }
};
