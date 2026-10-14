export default {
  template: "<div></div>",
  props: ["options"],
  async mounted() {
    this.alive = true;
    this.dirty = true;
    const { Plotly } = await import("nicegui-plotly");
    if (!this.alive) return;
    this.Plotly = Plotly;
    this.resizeObserver = new ResizeObserver(() => this.renderVisible());
    this.resizeObserver.observe(this.$el);
    this.renderVisible();
  },
  unmounted() {
    this.alive = false;
    this.resizeObserver?.disconnect();
    if (this.Plotly) this.Plotly.purge(this.$el);
  },
  methods: {
    update() {
      this.dirty = true;
      this.renderVisible();
    },
    async renderVisible() {
      const width = this.$el.clientWidth;
      if (!this.alive || !this.Plotly || this.rendering || width === 0 ||
          this.$el.getClientRects().length === 0) return;
      if (!this.dirty && width === this.lastWidth) return;
      this.dirty = false;
      this.rendering = true;
      this.lastWidth = width;
      try {
        // Own the responsive sizing. Plotly's deferred window resize can fire
        // after a tab is hidden and reject its promise for a non-displayed div.
        await this.Plotly.react(this.$el, this.options.data,
          { ...this.options.layout, width },
          { ...this.options.config, responsive: false });
      } finally {
        this.rendering = false;
        if (this.alive) this.renderVisible();
      }
    },
    async run_plot_method(name, ...args) {
      const parts = name.split(".");
      let owner = this.Plotly;
      for (const part of parts.slice(0, -1)) owner = owner[part];
      const result = await owner[parts.at(-1)](this.$el, ...args);
      return result === this.$el ? undefined : result;
    },
  },
};
