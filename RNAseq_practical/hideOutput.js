$(document).ready(function() {
  $output = $("pre").not(".r");
  $output.prepend("<div class='.showopt'>Show Output</div><br/>");
  $output.children("code").css({display: "none"});
            
  $(".showopt").click(function() {
    $btn = $(this);
    $chunk = $(this).parent().children("code");
    if($btn.html() === "Show Output") {
      $btn.html("Hide Output");
    } else {
      $btn.html("Show Output");
    }
    $chunk.slideToggle("fast", "swing");
  });
});
