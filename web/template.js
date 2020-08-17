window.Template = (function(){
    let templates = {};

    document.querySelectorAll('[type="text/template"]').forEach((e) => {
      templates[e.id] = e.innerHTML;
    });

    function replace(tpl, data){
        if (!tpl) return '';

        // "{:tag_name}" -> "tag_name"
        let trim = function(s){ return s.substr(2, s.length - 3); };

        // Variables
        // {:name}
        (tpl.match(/\{:([A-z0-9]+)*\}/g) || []).map(function(tag){
            let t = trim(tag);
            let r = new RegExp("{:" + t + "}", "g");
            tpl = tpl.replace(r, data[t] == undefined ? '' : data[t]);
        });

        // Includes
        // {:name use template}
        (tpl.match(/\{:[A-z0-9]+\ use [A-z-]+\}/g) || []).map(function(tag){
            let inc = trim(tag).split(' use ');
            let r = new RegExp(tag, "g");
            tpl = tpl.replace(r, Template(inc[1], data[inc[0]] || []));
        });

        return tpl;
    }

    return function(name, content){
        if (content.length != undefined) {
            return content.map(function(obj){
                return replace(templates[name], obj);
            }).join('');
        }
        return replace(templates[name], content);
    }
})();

/*

# SOURCE:
https://gist.github.com/latur/d67fde486b23312efbeac3d3e25582c2

# USAGE:
## Define templates (html):
<script type="text/template" id="user">
    <p>Hello, {:name}</p>
</script>
<script type="text/template" id="items">
    <p>Items: {:items use item-template}</p>
</script>
<script type="text/template" id="item-template">
    <b>Item: {:name}; <small>value: {:value}</small></b>
</script>
## Get html (JS):
window.Template = require('../components/template.js');;
> Template('user', {name: 'World'});
Result:
~~~
<p>Hello, World</p>
~~~
> Template('items', {
  title: 'My list',
  items: [
    {name: 'First', value: 23},
    {name: 'Second', value: 0}
  ]
});
Result:
~~~
<p>Items:
  <b>Item: First; <small>value: 23</small></b>
  <b>Item: Second; <small>value: 0</small></b>
</p>
~~~
*/
