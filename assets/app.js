import { firebaseConfig } from "./firebase-config.js";
import { initializeApp } from "https://www.gstatic.com/firebasejs/10.12.2/firebase-app.js";
import {
  getAuth, signInAnonymously, onAuthStateChanged
} from "https://www.gstatic.com/firebasejs/10.12.2/firebase-auth.js";
import {
  getFirestore, collection, addDoc, getDocs, query, orderBy, serverTimestamp
} from "https://www.gstatic.com/firebasejs/10.12.2/firebase-firestore.js";

// ---------------- Program data ----------------
const PROGRAM = [
  {
    id: "A",
    title: "Day A — Upper Back & Shoulders",
    sub: "Rear delts, mid-traps, rhomboids",
    exercises: [
      { name: "Band Pull-Apart (eccentric)", detail: "3 x 10" },
      { name: "DB Bent-Over Reverse Fly", detail: "3 x 8" },
      { name: "Wall Slides", detail: "3 x 8" },
      { name: "Band Face Pull", detail: "3 x 10" },
      { name: "Plank + Scap Squeeze", detail: "3 x 20-30s" },
      { name: "Doorway Chest Stretch", detail: "2 x 30s/side" },
    ],
  },
  {
    id: "B",
    title: "Day B — Posterior Chain & Core",
    sub: "Erectors, glutes, hamstrings",
    exercises: [
      { name: "Superman", detail: "3 x 8" },
      { name: "Bird Dog", detail: "3 x 8/side" },
      { name: "DB Romanian Deadlift", detail: "3 x 8" },
      { name: "Glute Bridge", detail: "3 x 10" },
      { name: "Band Row", detail: "3 x 10" },
      { name: "Kneeling Hip Flexor Stretch", detail: "2 x 30s/side" },
    ],
  },
  {
    id: "C",
    title: "Day C — Full Postural Integration",
    sub: "Everything working together",
    exercises: [
      { name: "Y-T-W Raises", detail: "3 x 6/letter" },
      { name: "DB Shrug", detail: "3 x 10" },
      { name: "Band Pull-Apart (high)", detail: "3 x 10" },
      { name: "Side Plank", detail: "2 x 15-20s/side" },
      { name: "Cat-Cow", detail: "2 x 8" },
      { name: "Child's Pose Reach", detail: "2 x 30s" },
    ],
  },
];

// ---------------- Firebase init ----------------
const app = initializeApp(firebaseConfig);
const auth = getAuth(app);
const db = getFirestore(app);
let sessions = [];

const authState = document.getElementById("authState");

onAuthStateChanged(auth, async (user) => {
  if (user) {
    authState.textContent = "connected";
    await loadSessions();
    render();
  }
});
signInAnonymously(auth).catch((err) => {
  authState.textContent = "connection failed";
  console.error(err);
});

// ---------------- Data ----------------
async function loadSessions() {
  const q = query(collection(db, "sessions"), orderBy("timestamp", "desc"));
  const snap = await getDocs(q);
  sessions = snap.docs.map((d) => {
    const data = d.data();
    return { ...data, date: data.timestamp?.toDate?.() ?? new Date() };
  });
}

async function logSession(dayId) {
  await addDoc(collection(db, "sessions"), {
    day: dayId,
    timestamp: serverTimestamp(),
  });
  await loadSessions();
  render();
}

// ---------------- Date helpers ----------------
function startOfWeek(d) {
  const date = new Date(d);
  const day = (date.getDay() + 6) % 7; // Monday = 0
  date.setDate(date.getDate() - day);
  date.setHours(0, 0, 0, 0);
  return date;
}

function sessionsInWeek(weekStart) {
  const weekEnd = new Date(weekStart);
  weekEnd.setDate(weekEnd.getDate() + 7);
  return sessions.filter((s) => s.date >= weekStart && s.date < weekEnd);
}

// ---------------- Rendering ----------------
const checkedState = {}; // { dayId: Set of exercise indices }
PROGRAM.forEach((d) => (checkedState[d.id] = new Set()));

function renderProgram() {
  const container = document.getElementById("program");
  container.innerHTML = "";

  const thisWeekStart = startOfWeek(new Date());
  const loggedToday = sessions.some((s) => {
    const t = new Date(s.date);
    const now = new Date();
    return t.toDateString() === now.toDateString();
  });

  PROGRAM.forEach((day) => {
    const card = document.createElement("div");
    card.className = "day-card";

    const allChecked = checkedState[day.id].size === day.exercises.length;

    card.innerHTML = `
      <div class="day-card__head">
        <h3>${day.title}</h3>
        <span class="day-tag">${day.exercises.length} MOVEMENTS</span>
      </div>
      <p class="day-card__sub">${day.sub}</p>
      <div class="exercise-list">
        ${day.exercises
          .map(
            (ex, i) => `
          <label class="exercise-row">
            <input type="checkbox" data-day="${day.id}" data-idx="${i}" ${
              checkedState[day.id].has(i) ? "checked" : ""
            } />
            <span class="exercise-name ${
              checkedState[day.id].has(i) ? "checked-text" : ""
            }">${ex.name}</span>
            <span class="exercise-detail">${ex.detail}</span>
          </label>`
          )
          .join("")}
      </div>
      <button class="log-btn" data-day="${day.id}" ${
        allChecked ? "" : "disabled"
      }>
        ${allChecked ? "Log " + day.id + " as complete" : "Check off all movements to log"}
      </button>
    `;

    container.appendChild(card);
  });

  container.querySelectorAll('input[type="checkbox"]').forEach((box) => {
    box.addEventListener("change", (e) => {
      const dayId = e.target.dataset.day;
      const idx = parseInt(e.target.dataset.idx, 10);
      if (e.target.checked) checkedState[dayId].add(idx);
      else checkedState[dayId].delete(idx);
      renderProgram();
    });
  });

  container.querySelectorAll(".log-btn").forEach((btn) => {
    btn.addEventListener("click", async (e) => {
      const dayId = e.target.dataset.day;
      btn.disabled = true;
      btn.textContent = "Logging…";
      await logSession(dayId);
      checkedState[dayId].clear();
    });
  });
}

function renderPlumb() {
  const thisWeekStart = startOfWeek(new Date());
  const count = sessionsInWeek(thisWeekStart).length;
  const target = 3;
  const deficit = Math.max(0, target - count);
  const angle = Math.min(deficit * 7, 18); // degrees, straightens as sessions increase

  document.getElementById("plumbString").style.transform = `rotate(${angle}deg)`;
  document.getElementById("plumbBob").style.transform = `rotate(${angle}deg)`;

  const label = document.getElementById("plumbLabel");
  if (count >= target) label.textContent = "plumb — on target";
  else if (count === 0) label.textContent = "off plumb — 0 this week";
  else label.textContent = `${target - count} session${target - count > 1 ? "s" : ""} to plumb`;
}

function renderStats() {
  const thisWeekStart = startOfWeek(new Date());
  document.getElementById("weekCount").textContent = sessionsInWeek(thisWeekStart).length;
  document.getElementById("totalCount").textContent = sessions.length;

  // streak = consecutive weeks (ending this week) with >= 3 sessions
  let streak = 0;
  let cursor = new Date(thisWeekStart);
  while (true) {
    const c = sessionsInWeek(cursor).length;
    if (c >= 3) {
      streak++;
      cursor.setDate(cursor.getDate() - 7);
    } else {
      // allow current week to still be "in progress" without breaking streak count of prior weeks
      if (cursor.getTime() === thisWeekStart.getTime()) {
        cursor.setDate(cursor.getDate() - 7);
        continue;
      }
      break;
    }
  }
  document.getElementById("streakCount").textContent = streak;
}

function renderHistory() {
  const grid = document.getElementById("historyGrid");
  grid.innerHTML = "";
  const thisWeekStart = startOfWeek(new Date());

  const weeks = [];
  for (let i = 9; i >= 0; i--) {
    const ws = new Date(thisWeekStart);
    ws.setDate(ws.getDate() - i * 7);
    weeks.push(ws);
  }

  weeks.forEach((ws) => {
    const count = sessionsInWeek(ws).length;
    const fill = count >= 3 ? "full" : count >= 1 ? "half" : "empty";
    const box = document.createElement("div");
    box.className = "week-box";
    box.dataset.fill = fill;
    box.title = `Week of ${ws.toLocaleDateString()}: ${count} session(s)`;
    grid.appendChild(box);
  });
}

function render() {
  renderProgram();
  renderPlumb();
  renderStats();
  renderHistory();
}

render();