$('.tissue_click').click(function(){
    var input = $( "#tissueFormOutput" );
    var add_text = this.innerHTML
    input.val( input.val() + add_text + " " );
 })


$('.tf_click').click(function(){
    var input = $( "#tfFormOutput" );
    var add_text = this.innerHTML
    input.val( input.val() + add_text + " " );
 })