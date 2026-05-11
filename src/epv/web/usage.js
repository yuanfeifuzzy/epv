
(() => {
  const MODAL_ID = 'helpModal';
  const TARGET_ID = 'helpModalBody';
  const README = 'README.md'; // same-origin path to your README

  // Optional: minimal styles for code blocks when injecting Markdown
  const style = document.createElement('style');
  style.textContent = `
    #${TARGET_ID} h1, #${TARGET_ID} h2, #${TARGET_ID} h3 { scroll-margin-top: 72px; }
    #${TARGET_ID} pre code { display:block; padding:1rem; background:#f6f8fa; border-radius:.5rem; overflow:auto; }
    #${TARGET_ID} code { background:#f6f8fa; padding:.15rem .35rem; border-radius:.25rem; }
  `;
  document.head.appendChild(style);

  async function loadReadmeInto(container) {
    container.innerHTML = '<div class="text-secondary">Loading README…</div>';
    try {
      // Cache-bust in dev so updates show immediately
      const res = await fetch(`${README}?t=${Date.now()}`, { credentials: 'same-origin' });
      if (!res.ok) throw new Error(`${res.status} ${res.statusText}`);
      const md = await res.text();

      // Render Markdown -> HTML
      const html = marked.parse(md, { mangle: false, headerIds: true });

      // Optional: sanitize if your README contains HTML (uncomment if you include DOMPurify)
      // const safe = DOMPurify.sanitize(html);
      // container.innerHTML = safe;

      container.innerHTML = html;

      // Optional: handle in-page hash links inside the modal
      container.querySelectorAll('a[href^="#"]').forEach(a => {
        a.addEventListener('click', (e) => {
          const id = a.getAttribute('href').slice(1);
          const target = container.querySelector(`#${CSS.escape(id)}`);
          if (target) {
            e.preventDefault();
            target.scrollIntoView({ behavior: 'smooth', block: 'start' });
          }
        });
      });
    } catch (err) {
      container.innerHTML = `
        <div class="alert alert-danger mb-0">
          Failed to load README.md: <code>${String(err.message || err)}</code>
        </div>`;
    }
  }

  document.addEventListener('DOMContentLoaded', () => {
    const modalEl = document.getElementById(MODAL_ID);
    const bodyEl  = document.getElementById(TARGET_ID);
    if (!modalEl || !bodyEl) return;

    // Load once on first open; set { once:false } if you want to reload every time
    modalEl.addEventListener('show.bs.modal', () => {
      if (bodyEl.dataset.loaded === '1') return;
      loadReadmeInto(bodyEl).then(() => { bodyEl.dataset.loaded = '1'; });
    }, { once: false });
  });
})();
