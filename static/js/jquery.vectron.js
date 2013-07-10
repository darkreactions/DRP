/*!
 * jQuery Vectron: jQuery plugin to load external SVG files into Raphael JS through declarative markup.
 * 
 * Copyright  2011 Room & Board
 * MIT Licensed
 */
(function($){
	
	$.vectron = {
		cache: {},
		defaults: {
			scale: 1
		}
	};
	
	$.fn.vectron = function(options) {
		var settings = $.extend({}, $.vectron.defaults, options);
		
		this.each(function(){
			if ( $(this).data('vectron') ){
				return;
			}
			
			$(this)
				.data('vectron', {
					src: $(this).attr('data-svg'),
					paper: Raphael(this, $(this).width(), $(this).height()),
					options: settings
				})
				.on('import.vectron', function(e){
					var $self = $(this),
						data = $self.data('vectron'),
						promise, svgText;
						
					if ( $.vectron.cache[data.src] ) {
						svgText = $.vectron.cache[data.src];		
						data.vectorItems = rappar(svgText).slice(0);
						$self.trigger('imported.vectron').trigger('render.vectron');
						return this;
					}
					
					promise = $.ajax({ url: data.src, type: "GET", dataType: "text" });
					promise.success(function(resp){
						if (!resp) return;
						
						svgText = resp.slice(resp.indexOf('<svg'));
						
						if (svgText) {
							var data = $self.data('vectron');
							$.vectron.cache[data.src] = svgText;
							data.vectorItems = rappar(svgText).slice(0);
							$self.trigger('imported.vectron').trigger('render.vectron');
						}
					});
					
					return this;
				})
				.on('render.vectron', function(e){
					var data = $(this).data('vectron');
					
					data.set = data.paper.add(data.vectorItems);
					$(this).trigger('afterrender.vectron').trigger('setscale.vectron');
				})
				.on('setscale.vectron', function(e) {
					var data = $(this).data('vectron');
					
					if (data.options.scale != 1) {
						data.set.transform("...S" + data.options.scale + "," + data.options.scale + ",0,0");
					}
					$(this).trigger('aftersetscale.vectron').trigger('setdimensions.vectron');
				})
				.on('setdimensions.vectron', function(e) {
					var data = $(this).data('vectron'),
						bbox = data.set.getBBox();
						
					data.paper.setSize(bbox.width + 3, bbox.height);
					$(this).width(bbox.width + 3).trigger('complete.vectron');
				});
				
		});
		
		return this.trigger('import.vectron'); // before returning, trigger import
	}
	
})(jQuery);